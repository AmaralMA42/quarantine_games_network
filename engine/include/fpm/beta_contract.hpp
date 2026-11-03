#pragma once
#include <cmath>
#include <stdexcept>
#include <string>
namespace fpm {
enum class BetaMixingRule { INTERPOLATED, LEGACY_PUBLISHED };
inline BetaMixingRule parse_beta_mixing_rule(const std::string &s){
 if(s=="INTERPOLATED") return BetaMixingRule::INTERPOLATED;
 if(s=="LEGACY_PUBLISHED") return BetaMixingRule::LEGACY_PUBLISHED;
 throw std::runtime_error("unknown beta_mixing_rule");
}
inline const char* beta_mixing_rule_name(BetaMixingRule r){return r==BetaMixingRule::INTERPOLATED?"INTERPOLATED":"LEGACY_PUBLISHED";}
struct ResolvedBetas { double beta_n,beta_a,beta_q,a; BetaMixingRule rule; bool claim_grade_corrected_model; };
inline ResolvedBetas resolve_betas(double bn,double bq,double a,BetaMixingRule rule,bool allow_legacy=false){
 if(!std::isfinite(bn)||!std::isfinite(bq)||!std::isfinite(a)) throw std::runtime_error("non-finite beta contract input");
 if(!(bq>0)) throw std::runtime_error("beta_Q must be positive");
 if(!(bn>bq)) throw std::runtime_error("beta_N must exceed beta_Q");
 double ba;
 if(rule==BetaMixingRule::INTERPOLATED){
  if(!(a>0&&a<1)) throw std::runtime_error("INTERPOLATED requires 0<a<1");
  ba=bq+a*(bn-bq);
  if(!(bn>ba&&ba>bq)) throw std::runtime_error("beta hierarchy violated");
  return {bn,ba,bq,a,rule,true};
 }
 ba=a*(bn+bq)/2.0;
 if(!(bn>ba&&ba>bq)&&!allow_legacy) throw std::runtime_error("LEGACY hierarchy violation requires override");
 return {bn,ba,bq,a,rule,false};
}
}
