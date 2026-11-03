#pragma once

#include <array>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

namespace fpm {

/** Minimal SHA-256 implementation used for provenance validation at the C++ boundary. */
class Sha256 {
public:
    Sha256() { reset(); }

    void reset() {
        state_ = {0x6a09e667u,0xbb67ae85u,0x3c6ef372u,0xa54ff53au,
                  0x510e527fu,0x9b05688cu,0x1f83d9abu,0x5be0cd19u};
        bit_length_ = 0; buffer_length_ = 0;
    }

    void update(const unsigned char* data, std::size_t length) {
        for (std::size_t i=0;i<length;++i) {
            buffer_[buffer_length_++] = data[i];
            if (buffer_length_ == 64) {
                transform(buffer_.data());
                bit_length_ += 512;
                buffer_length_ = 0;
            }
        }
    }

    void update(const std::string& text) {
        update(reinterpret_cast<const unsigned char*>(text.data()), text.size());
    }

    std::array<unsigned char,32> digest() {
        auto saved_state=state_; auto saved_buffer=buffer_;
        auto saved_bit_length=bit_length_; auto saved_buffer_length=buffer_length_;

        std::size_t i=buffer_length_;
        buffer_[i++]=0x80;
        if (i>56) {
            while (i<64) buffer_[i++]=0;
            transform(buffer_.data()); i=0;
        }
        while (i<56) buffer_[i++]=0;
        bit_length_ += static_cast<std::uint64_t>(buffer_length_)*8;
        for (int shift=56;shift>=0;shift-=8)
            buffer_[i++]=static_cast<unsigned char>((bit_length_>>shift)&0xffu);
        transform(buffer_.data());

        std::array<unsigned char,32> out{};
        for (std::size_t word=0;word<8;++word) {
            out[4*word]=static_cast<unsigned char>((state_[word]>>24)&0xffu);
            out[4*word+1]=static_cast<unsigned char>((state_[word]>>16)&0xffu);
            out[4*word+2]=static_cast<unsigned char>((state_[word]>>8)&0xffu);
            out[4*word+3]=static_cast<unsigned char>(state_[word]&0xffu);
        }
        state_=saved_state; buffer_=saved_buffer;
        bit_length_=saved_bit_length; buffer_length_=saved_buffer_length;
        return out;
    }

    std::string hexdigest() {
        const auto bytes=digest();
        std::ostringstream out; out<<std::hex<<std::setfill('0');
        for (unsigned char byte:bytes) out<<std::setw(2)<<static_cast<unsigned int>(byte);
        return out.str();
    }

private:
    static constexpr std::array<std::uint32_t,64> K_ = {
        0x428a2f98u,0x71374491u,0xb5c0fbcfu,0xe9b5dba5u,0x3956c25bu,0x59f111f1u,0x923f82a4u,0xab1c5ed5u,
        0xd807aa98u,0x12835b01u,0x243185beu,0x550c7dc3u,0x72be5d74u,0x80deb1feu,0x9bdc06a7u,0xc19bf174u,
        0xe49b69c1u,0xefbe4786u,0x0fc19dc6u,0x240ca1ccu,0x2de92c6fu,0x4a7484aau,0x5cb0a9dcu,0x76f988dau,
        0x983e5152u,0xa831c66du,0xb00327c8u,0xbf597fc7u,0xc6e00bf3u,0xd5a79147u,0x06ca6351u,0x14292967u,
        0x27b70a85u,0x2e1b2138u,0x4d2c6dfcu,0x53380d13u,0x650a7354u,0x766a0abbu,0x81c2c92eu,0x92722c85u,
        0xa2bfe8a1u,0xa81a664bu,0xc24b8b70u,0xc76c51a3u,0xd192e819u,0xd6990624u,0xf40e3585u,0x106aa070u,
        0x19a4c116u,0x1e376c08u,0x2748774cu,0x34b0bcb5u,0x391c0cb3u,0x4ed8aa4au,0x5b9cca4fu,0x682e6ff3u,
        0x748f82eeu,0x78a5636fu,0x84c87814u,0x8cc70208u,0x90befffau,0xa4506cebu,0xbef9a3f7u,0xc67178f2u};
    std::array<std::uint32_t,8> state_{};
    std::array<unsigned char,64> buffer_{};
    std::uint64_t bit_length_=0;
    std::size_t buffer_length_=0;

    static std::uint32_t rotr(std::uint32_t x, unsigned n) { return (x>>n)|(x<<(32-n)); }
    void transform(const unsigned char* block) {
        std::array<std::uint32_t,64> w{};
        for (int i=0;i<16;++i)
            w[i]=(std::uint32_t(block[4*i])<<24)|(std::uint32_t(block[4*i+1])<<16)|
                 (std::uint32_t(block[4*i+2])<<8)|std::uint32_t(block[4*i+3]);
        for (int i=16;i<64;++i) {
            const auto s0=rotr(w[i-15],7)^rotr(w[i-15],18)^(w[i-15]>>3);
            const auto s1=rotr(w[i-2],17)^rotr(w[i-2],19)^(w[i-2]>>10);
            w[i]=w[i-16]+s0+w[i-7]+s1;
        }
        auto a=state_[0],b=state_[1],c=state_[2],d=state_[3];
        auto e=state_[4],f=state_[5],g=state_[6],h=state_[7];
        for (int i=0;i<64;++i) {
            const auto S1=rotr(e,6)^rotr(e,11)^rotr(e,25);
            const auto ch=(e&f)^((~e)&g);
            const auto temp1=h+S1+ch+K_[i]+w[i];
            const auto S0=rotr(a,2)^rotr(a,13)^rotr(a,22);
            const auto maj=(a&b)^(a&c)^(b&c);
            const auto temp2=S0+maj;
            h=g; g=f; f=e; e=d+temp1; d=c; c=b; b=a; a=temp1+temp2;
        }
        state_[0]+=a; state_[1]+=b; state_[2]+=c; state_[3]+=d;
        state_[4]+=e; state_[5]+=f; state_[6]+=g; state_[7]+=h;
    }
};

/** Return lowercase SHA-256 for a byte string. */
inline std::string sha256_hex(const std::string& payload) {
    Sha256 sha; sha.update(payload); return sha.hexdigest();
}

/** Return lowercase SHA-256 for a file, failing closed on read errors. */
inline std::string sha256_file(const std::string& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("could not open file for SHA-256: " + path);
    Sha256 sha; std::array<char,1<<16> buffer{};
    while (input) {
        input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        const auto count=input.gcount();
        if (count>0) sha.update(reinterpret_cast<const unsigned char*>(buffer.data()), static_cast<std::size_t>(count));
    }
    if (!input.eof()) throw std::runtime_error("I/O error while hashing: " + path);
    return sha.hexdigest();
}

} // namespace fpm
