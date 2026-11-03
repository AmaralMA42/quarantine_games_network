#pragma once

#include <cctype>
#include <charconv>
#include <fstream>
#include <iterator>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "fpm/sha256.hpp"

namespace fpm {

struct FlatJsonValue {
    enum class Kind { String, Integer, Boolean, StringArray };
    Kind kind = Kind::String;
    std::string string_value;
    long long integer_value = 0;
    bool boolean_value = false;
    std::vector<std::string> array_value;
};

class FlatJsonObjectParser {
    const std::string &text_;
    size_t pos_ = 0;

    void whitespace() {
        while (pos_ < text_.size()
               && std::isspace(static_cast<unsigned char>(text_[pos_]))) ++pos_;
    }

    void require(char expected) {
        whitespace();
        if (pos_ >= text_.size() || text_[pos_] != expected)
            throw std::runtime_error("invalid periodic certificate JSON");
        ++pos_;
    }

    std::string string_value() {
        require('"');
        std::string out;
        while (pos_ < text_.size()) {
            const unsigned char ch = static_cast<unsigned char>(text_[pos_++]);
            if (ch == '"') return out;
            if (ch < 0x20) throw std::runtime_error("control byte in periodic certificate JSON string");
            if (ch != '\\') {
                if (ch > 0x7f) throw std::runtime_error("periodic certificate JSON must be ASCII");
                out.push_back(static_cast<char>(ch));
                continue;
            }
            if (pos_ >= text_.size()) throw std::runtime_error("truncated periodic certificate escape");
            const char escaped = text_[pos_++];
            switch (escaped) {
                case '"': out.push_back('"'); break;
                case '\\': out.push_back('\\'); break;
                case '/': out.push_back('/'); break;
                case 'b': out.push_back('\b'); break;
                case 'f': out.push_back('\f'); break;
                case 'n': out.push_back('\n'); break;
                case 'r': out.push_back('\r'); break;
                case 't': out.push_back('\t'); break;
                case 'u': {
                    if (pos_ + 4 > text_.size()) throw std::runtime_error("truncated periodic certificate unicode escape");
                    unsigned value = 0;
                    for (int i = 0; i < 4; ++i) {
                        const char digit = text_[pos_++];
                        value <<= 4;
                        if (digit >= '0' && digit <= '9') value += unsigned(digit - '0');
                        else if (digit >= 'a' && digit <= 'f') value += unsigned(digit - 'a' + 10);
                        else if (digit >= 'A' && digit <= 'F') value += unsigned(digit - 'A' + 10);
                        else throw std::runtime_error("invalid periodic certificate unicode escape");
                    }
                    if (value > 0x7f) throw std::runtime_error("periodic certificate JSON must be ASCII");
                    out.push_back(static_cast<char>(value));
                    break;
                }
                default: throw std::runtime_error("invalid periodic certificate escape");
            }
        }
        throw std::runtime_error("unterminated periodic certificate JSON string");
    }

    FlatJsonValue value() {
        whitespace();
        if (pos_ >= text_.size()) throw std::runtime_error("missing periodic certificate JSON value");
        FlatJsonValue out;
        if (text_[pos_] == '"') {
            out.kind = FlatJsonValue::Kind::String;
            out.string_value = string_value();
            return out;
        }
        if (text_[pos_] == '[') {
            out.kind = FlatJsonValue::Kind::StringArray;
            ++pos_; whitespace();
            if (pos_ < text_.size() && text_[pos_] == ']') { ++pos_; return out; }
            while (true) {
                out.array_value.push_back(string_value());
                whitespace();
                if (pos_ < text_.size() && text_[pos_] == ']') { ++pos_; return out; }
                require(',');
            }
        }
        if (text_.compare(pos_, 4, "true") == 0) {
            pos_ += 4; out.kind = FlatJsonValue::Kind::Boolean; out.boolean_value = true; return out;
        }
        if (text_.compare(pos_, 5, "false") == 0) {
            pos_ += 5; out.kind = FlatJsonValue::Kind::Boolean; out.boolean_value = false; return out;
        }
        const size_t begin = pos_;
        if (text_[pos_] == '-') ++pos_;
        while (pos_ < text_.size() && std::isdigit(static_cast<unsigned char>(text_[pos_]))) ++pos_;
        if (pos_ == begin || (text_[begin] == '-' && pos_ == begin + 1))
            throw std::runtime_error("unsupported periodic certificate JSON value");
        out.kind = FlatJsonValue::Kind::Integer;
        out.string_value = text_.substr(begin, pos_ - begin);
        return out;
    }

public:
    explicit FlatJsonObjectParser(const std::string &text) : text_(text) {}

    std::map<std::string, FlatJsonValue> parse() {
        std::map<std::string, FlatJsonValue> result;
        require('{'); whitespace();
        if (pos_ < text_.size() && text_[pos_] == '}') ++pos_;
        else {
            while (true) {
                const std::string key = string_value();
                require(':');
                if (!result.emplace(key, value()).second)
                    throw std::runtime_error("duplicate periodic certificate key: " + key);
                whitespace();
                if (pos_ < text_.size() && text_[pos_] == '}') { ++pos_; break; }
                require(',');
            }
        }
        whitespace();
        if (pos_ != text_.size()) throw std::runtime_error("trailing periodic certificate JSON data");
        return result;
    }
};

inline bool is_spatial_periodic_topology(const std::string &topology) {
    return topology == "honeycomb" || topology == "truncated_square_488"
        || topology == "4.8.8" || topology == "square"
        || topology == "kagome" || topology == "triangular";
}

struct PeriodicCertificateEvidence {
    std::string sha256;
    int lx = 0;
    int ly = 0;
};

inline PeriodicCertificateEvidence validate_periodic_certificate(
    const std::string &certificate_path,
    const std::string &expected_certificate_sha256,
    const std::string &topology,
    const std::string &graph_sha256
) {
    std::ifstream input(certificate_path, std::ios::binary);
    if (!input) throw std::runtime_error("periodic topology certificate missing");
    const std::string bytes(
        (std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>()
    );
    const std::string actual_sha256 = sha256_hex(bytes);
    if (actual_sha256 != expected_certificate_sha256)
        throw std::runtime_error("periodic certificate SHA-256 mismatch");
    const auto object = FlatJsonObjectParser(bytes).parse();
    auto string_field = [&](const std::string &key) -> const std::string & {
        const auto found = object.find(key);
        if (found == object.end() || found->second.kind != FlatJsonValue::Kind::String)
            throw std::runtime_error("periodic certificate string field missing: " + key);
        return found->second.string_value;
    };
    auto integer_field = [&](const std::string &key) -> long long {
        const auto found = object.find(key);
        if (found == object.end() || found->second.kind != FlatJsonValue::Kind::Integer)
            throw std::runtime_error("periodic certificate integer field missing: " + key);
        long long parsed = 0;
        const auto &raw = found->second.string_value;
        const auto result = std::from_chars(raw.data(), raw.data() + raw.size(), parsed);
        if (result.ec != std::errc{} || result.ptr != raw.data() + raw.size())
            throw std::runtime_error("periodic certificate integer field invalid: " + key);
        return parsed;
    };
    auto boolean_field = [&](const std::string &key) -> bool {
        const auto found = object.find(key);
        if (found == object.end() || found->second.kind != FlatJsonValue::Kind::Boolean)
            throw std::runtime_error("periodic certificate boolean field missing: " + key);
        return found->second.boolean_value;
    };
    auto array_field = [&](const std::string &key) -> const std::vector<std::string> & {
        const auto found = object.find(key);
        if (found == object.end() || found->second.kind != FlatJsonValue::Kind::StringArray)
            throw std::runtime_error("periodic certificate string-array field missing: " + key);
        return found->second.array_value;
    };
    if (string_field("schema") != "fpm.periodic_graph.v1")
        throw std::runtime_error("periodic certificate schema mismatch");
    if (string_field("topology") != topology)
        throw std::runtime_error("periodic certificate topology mismatch");
    if (string_field("graph_sha256") != graph_sha256)
        throw std::runtime_error("periodic certificate graph SHA-256 mismatch");
    if (string_field("boundary_condition") != "periodic_torus")
        throw std::runtime_error("periodic certificate boundary mismatch");
    if (string_field("construction").empty())
        throw std::runtime_error("periodic certificate construction provenance missing");
    if (array_field("periodic_axes") != std::vector<std::string>{"x", "y"})
        throw std::runtime_error("periodic certificate must close x and y axes");
    if (!boolean_field("periodic_wrap_verified"))
        throw std::runtime_error("periodic wrap-around is not verified");
    const long long lx = integer_field("lx");
    const long long ly = integer_field("ly");
    if (lx <= 0 || ly <= 0 || lx > 1000000000LL || ly > 1000000000LL)
        throw std::runtime_error("periodic certificate dimensions invalid");
    if ((topology == "truncated_square_488" || topology == "4.8.8")
        && (lx < 4 || ly < 4 || lx % 2 || ly % 2))
        throw std::runtime_error("claim-grade 4.8.8 requires even periodic dimensions >=4");
    return {actual_sha256, static_cast<int>(lx), static_cast<int>(ly)};
}

} // namespace fpm
