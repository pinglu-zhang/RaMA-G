#include "detail.hpp"
#include "centromere/refinement/block_dp.hpp"
#include "centromere/stitch/kernel.hpp"
#include "centromere/stitch/fixed_gaps.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <map>
#include <new>
#include <set>
#include <sstream>
#include <tuple>
#include <unordered_map>
#if RAMAG_HAVE_CENTROMERE
#include "anchor.h"
#endif

namespace ramag {
using namespace centromere_refinement_detail;
namespace {
#if RAMAG_HAVE_CENTROMERE
// One extra effective rarity tier supplies a competing repeat placement. The
// cap of 100 remains a search bound, not a reason to enumerate all copy pairs.
struct JointEnoughLevels {};
struct JointSeed { Position r{},q{}; Length n{}; std::uint64_t rc{},qc{}; };
std::vector<JointSeed> JointSeeds(const SequenceRecord &r,const SequenceRecord &q,
    const CentromereEndpointReview &v,const CentromereStitchOptions &o) {
  const auto rs=r.bases.substr(v.reference_begin,v.reference_end-v.reference_begin);
  const auto qs=q.bases.substr(v.query_begin,v.query_end-v.query_begin);
  if(o.gap_memory_limit_bytes < 16ULL*1024*1024)
    throw detail::StitchGapResourceLimit("gap_memory_budget");
  std::vector<rama_source::SequenceInfo> data{
    {rs,"joint-reference",static_cast<uint_t>(rs.size())},
    {qs,"joint-query",static_cast<uint_t>(qs.size())}};
  rama_source::AnchorFinder finder(data,"",0,false,false,100,[&]{
    if(o.interruption_callback)o.interruption_callback("joint-local-candidates");},true);
  std::map<std::tuple<Position,Position,Length>,JointSeed> unique;
  std::set<std::uint64_t> levels;
  std::uint64_t visited=0;
  try {
    finder.visitRegionCandidates({0,static_cast<uint_t>(rs.size()),0,static_cast<uint_t>(qs.size())},true,
      [&](const rama_source::RareCandidate &c){
        if(o.interruption_callback)o.interruption_callback("joint-local-candidate");
        if(!levels.contains(c.level) && levels.size()==2)throw JointEnoughLevels{};
        levels.insert(c.level);
        if(++visited>100000)throw detail::StitchGapResourceLimit("candidate_budget");
        const auto &a=c.pair;
        if(a.match_length<64)return;
        if(a.second_pos<rs.size()+1)throw AlignmentError("joint candidate offset");
        const Position qp=a.second_pos-rs.size()-1;
        if(a.first_pos>rs.size() || a.match_length>rs.size()-a.first_pos ||
           qp>qs.size() || a.match_length>qs.size()-qp || !c.reference_count || !c.query_count)
          throw AlignmentError("joint candidate bounds/count");
        unique.emplace(std::tuple{a.first_pos,qp,a.match_length},
          JointSeed{a.first_pos,qp,a.match_length,c.reference_count,c.query_count});
        if(unique.size()>3000)throw detail::StitchGapResourceLimit("candidate_budget");
      });
  } catch(const JointEnoughLevels &) {}
  std::vector<JointSeed> a;
  for(const auto &[key,x]:unique){(void)key;a.push_back(x);}
  auto gap=[](Length rlen,Length qlen){
    const auto n=static_cast<std::int64_t>(rlen>qlen?rlen-qlen:qlen-rlen);
    return n?static_cast<double>(std::min({60+30*n,800+5*n,2500+n})):0.;
  };
  std::vector<double> score(a.size());
  std::vector<std::size_t> parent(a.size(),a.size());
  std::vector<Length> trim(a.size());
  double best=-std::numeric_limits<double>::infinity();std::size_t last=a.size();
  for(std::size_t i=0;i<a.size();++i){
    if(o.interruption_callback)o.interruption_callback("joint-seed-chain");
    const auto x=a[i];
    const double weight=20.*static_cast<double>(x.n)/
      std::sqrt(static_cast<double>(x.rc)*static_cast<double>(x.qc));
    score[i]=weight-gap(x.r,x.q);
    for(std::size_t j=0;j<i;++j){
      const auto y=a[j];
      const Length t=std::max(y.r+y.n>x.r?y.r+y.n-x.r:0,
                              y.q+y.n>x.q?y.q+y.n-x.q:0);
      if(t>=x.n || y.r+y.n>x.r+t || y.q+y.n>x.q+t)continue;
      const double value=score[j]+weight*static_cast<double>(x.n-t)/static_cast<double>(x.n)
        -gap(x.r+t-y.r-y.n,x.q+t-y.q-y.n);
      if(value>score[i]){score[i]=value;parent[i]=j;trim[i]=t;}
    }
    const auto value=score[i]-gap(rs.size()-x.r-x.n,qs.size()-x.q-x.n);
    if(value>best){best=value;last=i;}
  }
  std::vector<JointSeed> chain;
  for(auto i=last;i<a.size();i=parent[i]){
    auto x=a[i];x.r+=trim[i];x.q+=trim[i];x.n-=trim[i];
    // The maximal match discovers a locus; it must not freeze every repeat
    // base out to its endpoints. Preserve its occurrence class using minimal
    // prefix/suffix contexts, leaving ambiguous extensions to gap alignment.
    const auto counts=finder.countFullOccurrences(static_cast<uint_t>(x.r),static_cast<uint_t>(x.n));
    auto same=[&](Position begin,Length n){return finder.countFullOccurrences(static_cast<uint_t>(begin),static_cast<uint_t>(n))==counts;};
    Length lo=std::min<Length>(32,x.n),hi=x.n;
    while(lo<hi){const auto mid=lo+(hi-lo)/2;if(same(x.r,mid))hi=mid;else lo=mid+1;}
    const auto prefix=lo;
    lo=0;hi=prefix-std::min<Length>(32,prefix);
    while(lo<hi){const auto mid=lo+(hi-lo+1)/2;if(same(x.r+mid,prefix-mid))lo=mid;else hi=mid-1;}
    const auto begin=lo;
    lo=std::min<Length>(32,x.n);hi=x.n;
    while(lo<hi){const auto mid=lo+(hi-lo)/2;if(same(x.r+x.n-mid,mid))hi=mid;else lo=mid+1;}
    const auto suffix=lo,right_begin=x.n-suffix;
    lo=std::min<Length>(32,suffix);hi=suffix;
    while(lo<hi){const auto mid=lo+(hi-lo)/2;if(same(x.r+right_begin,mid))hi=mid;else lo=mid+1;}
    const auto context_first=std::min(begin,right_begin),context_last=std::max(prefix,right_begin+lo);
    chain.push_back({x.r+context_first+v.reference_begin,x.q+context_first+v.query_begin,context_last-context_first,x.rc,x.qc});

  }
  std::reverse(chain.begin(),chain.end());return chain;
}
// A clean set of old fragments supplies a conservative veto when a new
// path both weakens score per exact base and introduces substitution-like
// breaks. Already discordant sparse fragments do not calibrate that veto:
// their observed break density is not comparable to newly recovered coverage.
bool JointRetainsCopyRoot(const A&a,const CentromereEndpointReview&v){
  Position rp=v.reference_begin,qp=v.query_begin;Length retained=0;
  for(auto op:v.new_cigar){
    if(op.operation=='='){
      const auto lo=std::max(qp,a.query_begin),hi=std::min(qp+op.length,a.query_begin+a.length);
      if(lo<hi && rp+lo-qp==a.reference_begin+lo-a.query_begin)retained+=hi-lo;
    }
    if(op.operation!='I')rp+=op.length;
    if(op.operation!='D')qp+=op.length;
  }
  return retained==a.length;
}
// Union intervals on the same diagonal, or reject competing physical pairs.
bool JointMergeCollinear(std::vector<JointSeed>&seeds){
  std::sort(seeds.begin(),seeds.end(),[](const auto&a,const auto&b){
    return std::tuple{a.r,a.q}<std::tuple{b.r,b.q};
  });
  std::size_t count=0;
  for(const auto seed:seeds){
    if(count){
      auto& tail=seeds[count-1];
      if(seed.r<tail.r+tail.n || seed.q<tail.q+tail.n){
        if(seed.q<tail.q || seed.r-tail.r!=seed.q-tail.q)return false;
        tail.n=std::max(tail.n,seed.r-tail.r+seed.n);
        continue;
      }
    }
    seeds[count++]=seed;
  }
  seeds.resize(count);
  return true;
}

bool JointQualityDominated(const CrossFragmentSummary&old,std::int64_t score,Length exact,Length short_breaks){
  if(!old.exact||!exact)return false;
  // 65,536-base windows and disjoint primary query paths bound these integer
  // products; comparisons do not use rounded identities or floating epsilons.
  const auto before=static_cast<std::int64_t>(old.exact),after=static_cast<std::int64_t>(exact);
  return !old.short_disagreements && short_breaks && score*before<old.score*after;
}
bool JointUnique(const A &a,const SequenceRecord &r,const SequenceRecord &q){
  if(a.length<64)return false;
  const auto word=std::string_view(r.bases).substr(a.reference_begin,a.length);
  const auto &query=q.bases;
  const auto first_r=r.bases.find(word),first_q=query.find(word);
  return first_r==a.reference_begin && first_q==a.query_begin &&
    r.bases.find(word,first_r+1)==std::string::npos &&
    query.find(word,first_q+1)==std::string::npos;
}
// A long word can be unique only because it bridges two incompatible repeat
// ends. Require independently locating, disjoint variable-length prefix/suffix
// contexts, rather than promoting whole-word uniqueness to two witnesses.
bool JointVariableFlank(const A&a,const SequenceRecord&r,const SequenceRecord&q){
  if(!JointUnique(a,r,q))return false;
  const auto word=std::string_view(r.bases).substr(a.reference_begin,a.length);
  auto unique=[&](std::string_view text){
    const auto ri=r.bases.find(text),qi=q.bases.find(text);
    return ri!=std::string::npos&&qi!=std::string::npos&&r.bases.find(text,ri+1)==std::string::npos&&q.bases.find(text,qi+1)==std::string::npos;
  };
  auto minimum=[&](bool suffix){
    Length lo=32,hi=a.length;
    while(lo<hi){const auto mid=lo+(hi-lo)/2;
      if(unique(suffix?word.substr(a.length-mid):word.substr(0,mid)))hi=mid;else lo=mid+1;
    }
    return lo;
  };
  const auto prefix=minimum(false),suffix=minimum(true);
  return prefix<=a.length-suffix;
}
// If the candidate retained part of an independently positioned old flank,
// restore its ends only when the new trace leaves both coordinates free.
bool JointRestoreFlanks(CentromereEndpointReview&v,const std::vector<A>&anchors,
    const SequenceRecord&r,const SequenceRecord&q){
  std::vector<JointSeed> pairs;
  Position rp=v.reference_begin,qp=v.query_begin;
  for(auto op:v.new_cigar){
    if(op.operation=='=')pairs.push_back({rp,qp,op.length,0,0});
    if(op.operation!='I')rp+=op.length;
    if(op.operation!='D')qp+=op.length;
  }
  const auto original_size=pairs.size();
  for(const auto&a:anchors){
    if(a.strand!=Strand::Forward || a.query_begin<v.query_begin || a.query_begin+a.length>v.query_end ||
       a.reference_begin<v.reference_begin || a.reference_begin+a.length>v.reference_end ||
       !JointVariableFlank(a,r,q))continue;
    if(JointRetainsCopyRoot(a,v))continue;
    const bool retained=std::any_of(pairs.begin(),pairs.begin()+static_cast<std::ptrdiff_t>(original_size),
      [&](const auto&seed){
        const auto lo=std::max(seed.q,a.query_begin),hi=std::min(seed.q+seed.n,a.query_begin+a.length);
        return lo<hi && seed.r+lo-seed.q==a.reference_begin+lo-a.query_begin;
      });
    if(!retained)return false;
    pairs.push_back({a.reference_begin,a.query_begin,a.length,1,1});
  }
  if(pairs.size()==original_size)return true;
  if(!JointMergeCollinear(pairs))return false;
  v.new_cigar.clear();rp=v.reference_begin;qp=v.query_begin;
  for(const auto&seed:pairs){
    Emit(v.new_cigar,'D',seed.r-rp);Emit(v.new_cigar,'I',seed.q-qp);
    Emit(v.new_cigar,'=',seed.n);
    rp=seed.r+seed.n;qp=seed.q+seed.n;
  }
  Emit(v.new_cigar,'D',v.reference_end-rp);Emit(v.new_cigar,'I',v.query_end-qp);
  return true;
}

// A same-score sliding gap leaves its physical repeat copy unresolved. Short
// homopolymer endpoints retain the normal DP tie rule; a whole independent
// context (64 bases, as used by flank evidence) cannot be assigned by that tie.
bool JointMovableGap(const SequenceRecord&r,const SequenceRecord&q,const CentromereEndpointReview&v){
  Position rp=v.reference_begin,qp=v.query_begin;
  for(std::size_t k=0;k<v.new_cigar.size();++k){
    const auto op=v.new_cigar[k];
    if(op.operation=='I'||op.operation=='D'){
      const auto &sequence=op.operation=='I'?q.bases:r.bases;
      const Position pos=op.operation=='I'?qp:rp;
      const Length left=k && v.new_cigar[k-1].operation=='='?v.new_cigar[k-1].length:0;
      const Length right=k+1<v.new_cigar.size() && v.new_cigar[k+1].operation=='='?v.new_cigar[k+1].length:0;
      Length a=0,b=0;
      while(a<left && pos>a && sequence[pos-a-1]==sequence[pos+op.length-a-1])++a;
      while(b<right && pos+op.length+b<sequence.size() && sequence[pos+b]==sequence[pos+op.length+b])++b;
      if(a+b>=64)return true;
    }
    if(op.operation!='I')rp+=op.length;
    if(op.operation!='D')qp+=op.length;
  }
  return false;
}
void JointFilterCoincidences(CentromereEndpointReview&v,const CentromereStitchResult&old,const SequenceRecord&r,const SequenceRecord&q,const std::vector<Range>&confirmed){
  std::vector<Range> uncertain;
  Position arp=v.reference_begin,aqp=v.query_begin;
  for(std::size_t k=0;k<v.new_cigar.size();++k){
    const auto op=v.new_cigar[k];
    if(op.operation=='I'||op.operation=='D'){
      const auto &sequence=op.operation=='I'?q.bases:r.bases;
      const Position position=op.operation=='I'?aqp:arp;
      const auto left=k && v.new_cigar[k-1].operation=='='?v.new_cigar[k-1].length:0;
      const auto right=k+1<v.new_cigar.size()&&v.new_cigar[k+1].operation=='='?v.new_cigar[k+1].length:0;
      Length before=0,after=0;
      while(before+1<left&&position>before&&sequence[position-before-1]==sequence[position+op.length-before-1])++before;
      while(after+1<right&&position+op.length+after<sequence.size()&&sequence[position+after]==sequence[position+op.length+after])++after;
      if(before||after)uncertain.emplace_back(aqp-before,aqp+(op.operation=='I'?op.length:0)+after);
    }
    if(op.operation!='I')arp+=op.length;
    if(op.operation!='D')aqp+=op.length;
  }
  uncertain=Union(std::move(uncertain));
  std::vector<CigarOp> filtered;Position rp=v.reference_begin,qp=v.query_begin;
  for(auto op:v.new_cigar){
    if(op.operation!='=')Emit(filtered,op.operation,op.length);
    else{
      std::vector<Range> retained;
      for(const auto&segment:old.segments){
        const auto&a=segment.alignment;if(segment.path_role!="primary"||a.strand!=Strand::Forward)continue;
        Position rr=a.reference_begin,qq=a.query_begin;
        for(auto previous:a.cigar){
          if(previous.operation=='='){
            const auto lo=std::max(qp,qq),hi=std::min(qp+op.length,qq+previous.length);
            if(lo<hi&&rp+lo-qp==rr+lo-qq)retained.emplace_back(lo,hi);
          }
          if(previous.operation!='I')rr+=previous.length;
          if(previous.operation!='D')qq+=previous.length;
        }
      }
      Position end=qp;
      auto unpaired=[&](Length n){Emit(filtered,'D',n);Emit(filtered,'I',n);};
      auto confident=[&](Length n){if(n>=32)Emit(filtered,'=',n);else unpaired(n);};
      auto emit_new=[&](Position a,Position b){
        auto allowed=[&](Position begin,Position finish){
          for(auto [lo,hi]:uncertain){
            const auto first=std::max(begin,lo),last=std::min(finish,hi);
            if(first>=last)continue;
            confident(first-begin);unpaired(last-first);begin=last;
          }
          confident(finish-begin);
        };
        for(auto [lo,hi]:confirmed){
          const auto first=std::max(a,lo),last=std::min(b,hi);
          if(first>=last)continue;
          unpaired(first-a);allowed(first,last);a=last;
        }
        unpaired(b-a);
      };
      for(auto [a,b]:Union(std::move(retained))){emit_new(end,a);Emit(filtered,'=',b-a);end=b;}
      emit_new(end,qp+op.length);
    }
    if(op.operation!='I')rp+=op.length;
    if(op.operation!='D')qp+=op.length;
  }
  // Unpaired intervals carry no coordinate relation internally. Normalize
  // their total consumption so swapping inputs cannot change a gap-open cost
  // merely through the arbitrary order of newly rejected coincidences.
  std::vector<CigarOp> canonical;Length deleted=0,inserted=0;
  auto flush=[&]{Emit(canonical,'D',deleted);Emit(canonical,'I',inserted);deleted=inserted=0;};
  for(auto op:filtered){
    if(op.operation=='='){flush();Emit(canonical,'=',op.length);}
    else if(op.operation=='D')deleted+=op.length;
    else if(op.operation=='I')inserted+=op.length;
    else throw AlignmentError("joint confidence filter contains substitution");
  }
  flush();v.new_cigar=std::move(canonical);
}
bool JointLongGapWithoutContext(const SequenceRecord&r,const SequenceRecord&q,
    const CentromereStitchResult&old,const CentromereEndpointReview&v){
  Length unpaired=0;for(auto op:v.new_cigar)if(op.operation=='I'||op.operation=='D')unpaired+=op.length;
  if(unpaired<64)return false;
  Position rp=v.reference_begin,qp=v.query_begin;
  for(auto op:v.new_cigar){
    if(op.operation=='=' && op.length>=64){
      std::vector<Range> known;
      for(const auto&segment:old.segments){
        const auto&alignment=segment.alignment;
        if(segment.path_role!="primary"||alignment.strand!=Strand::Forward)continue;
        Position rr=alignment.reference_begin,qq=alignment.query_begin;
        for(auto previous:alignment.cigar){
          if(previous.operation=='='){
            const auto first=std::max(qp,qq),last=std::min(qp+op.length,qq+previous.length);
            if(first<last&&rp+first-qp==rr+first-qq)known.emplace_back(first,last);
          }
          if(previous.operation!='I')rr+=previous.length;
          if(previous.operation!='D')qq+=previous.length;
        }
      }
      auto new_context=[&](Position first,Position last){
        if(last-first<64)return false;
        const auto word=std::string_view(r.bases).substr(rp+first-qp,last-first);
        const auto rfirst=r.bases.find(word),qfirst=q.bases.find(word);
        return r.bases.find(word,rfirst+1)==std::string::npos&&q.bases.find(word,qfirst+1)==std::string::npos;
      };
      Position last=qp;
      for(auto[first,finish]:Union(std::move(known))){if(new_context(last,first))return false;last=finish;}
      if(new_context(last,qp+op.length))return false;
    }
    if(op.operation!='I')rp+=op.length;
    if(op.operation!='D')qp+=op.length;
  }
  return true;
}


// Keep all existing exact physical pairings while reopening their intervening
// gaps. If they are not a single ordered path this route is inapplicable: it
// must not silently choose which old copy or direction to discard.
std::optional<std::vector<JointSeed>> JointConstrainPairs(
    const CentromereStitchResult&old,const CentromereEndpointReview&v,
    const std::vector<JointSeed>&seeds){
  std::vector<JointSeed> fixed;
  for(const auto&segment:old.segments){
    const auto&a=segment.alignment;
    if(segment.path_role!="primary"||a.strand!=Strand::Forward)continue;
    Position rp=a.reference_begin,qp=a.query_begin;
    for(auto op:a.cigar){
      if(op.operation=='='||op.operation=='X'){
        const auto lo=std::max(v.reference_begin>rp?v.reference_begin-rp:0,v.query_begin>qp?v.query_begin-qp:0);
        const auto hi=std::min({op.length,v.reference_end>rp?v.reference_end-rp:0,v.query_end>qp?v.query_end-qp:0});
        if(lo<hi){if(op.operation!='=')return std::nullopt;fixed.push_back({rp+lo,qp+lo,hi-lo,0,0});}
      }
      if(op.operation!='I')rp+=op.length;
      if(op.operation!='D')qp+=op.length;
    }
  }
  if(fixed.empty())return std::nullopt;
  std::sort(fixed.begin(),fixed.end(),[](const auto&a,const auto&b){return std::tuple{a.q,a.r}<std::tuple{b.q,b.r};});
  std::vector<JointSeed> unique;
  for(const auto&a:fixed){
    if(!unique.empty()){
      auto&b=unique.back();
      if(a.q<b.q+b.n||a.r<b.r+b.n){
        if(a.q<b.q||a.r<b.r||a.q-b.q!=a.r-b.r)return std::nullopt;
        b.n=std::max(b.n,a.q+a.n-b.q);continue;
      }
    }
    unique.push_back(a);
  }
  std::vector<JointSeed> result;Position rb=v.reference_begin,qb=v.query_begin;
  for(std::size_t k=0;k<=unique.size();++k){
    const Position re=k<unique.size()?unique[k].r:v.reference_end;
    const Position qe=k<unique.size()?unique[k].q:v.query_end;
    for(const auto&a:seeds){
      const auto lo=std::max(rb>a.r?rb-a.r:0,qb>a.q?qb-a.q:0);
      const auto hi=std::min({a.n,re>a.r?re-a.r:0,qe>a.q?qe-a.q:0});
      if(lo<hi&&hi-lo>=32)result.push_back({a.r+lo,a.q+lo,hi-lo,a.rc,a.qc});
    }
    if(k<unique.size()){result.push_back(unique[k]);rb=unique[k].r+unique[k].n;qb=unique[k].q+unique[k].n;}
  }
  return result;
}
bool JointRetainsExactPairs(const CentromereStitchResult&old,const CentromereEndpointReview&v){
  const auto fixed=JointConstrainPairs(old,v,{});if(!fixed)return false;
  for(const auto&a:*fixed){
    Position rp=v.reference_begin,qp=v.query_begin;Length n=0;
    for(auto op:v.new_cigar){
      if(op.operation=='='){
        const auto lo=std::max(qp,a.q),hi=std::min(qp+op.length,a.q+a.n);
        if(lo<hi&&rp+lo-qp==a.r+lo-a.q)n+=hi-lo;
      }
      if(op.operation!='I')rp+=op.length;
      if(op.operation!='D')qp+=op.length;
    }
    if(n!=a.n)return false;
  }
  return true;
}


// A newly localized component with substitution-like opposing short gaps is
// not established by its exact islands alone. Withhold its new pairings while
// preserving every old physical pair; other independently bounded components
// of the same retained path can still be recovered.
void JointMaskDiscordantComponents(CentromereEndpointReview&v,
    const CentromereStitchResult&old,const CentromereDirectedAnchorResult&input){
  std::vector<Range> uncertain;Position rp=v.reference_begin,qp=v.query_begin;
  std::size_t discordances=0;bool terminal_disagreement=false;
  for(std::size_t i=0;i<v.new_cigar.size();++i){
    const auto a=v.new_cigar[i];
    if(i+1<v.new_cigar.size()){
      const auto b=v.new_cigar[i+1];
      if(((a.operation=='I'&&b.operation=='D')||(a.operation=='D'&&b.operation=='I'))&&a.length<64&&b.length<64){
        const auto consumed=a.operation=='I'?a.length:b.length;
        if(std::any_of(input.unknown.begin(),input.unknown.end(),[&](const auto&x){return Touch(qp,qp+consumed,x.query_begin,x.query_end);})){
          ++discordances;
          // The fixed exterior match does not establish a new endpoint across
          // a substitution-like terminal corner. Keep it unresolved until an
          // independent endpoint explanation is available on both sides.
          terminal_disagreement|=i==0||i+2==v.new_cigar.size()||
            std::any_of(input.unknown.begin(),input.unknown.end(),[&](const auto&x){
              return Touch(qp,qp+consumed,x.query_begin,x.query_end)&&(qp==x.query_begin||qp+consumed==x.query_end);
            });
        }
      }
    }
    if(a.operation!='I')rp+=a.length;
    if(a.operation!='D')qp+=a.length;
  }
  // One local indel corner can be ordinary endpoint normalization. Multiple
  // separate substitution-like corners contradict localization of the new
  // component; a short old island cannot split that contradiction away.
  if(discordances<2&&!terminal_disagreement)return;
  uncertain.emplace_back(v.query_begin,v.query_end);
  const auto fixed=JointConstrainPairs(old,v,{});
  if(!fixed)return;
  std::vector<CigarOp> filtered;rp=v.reference_begin;qp=v.query_begin;
  for(auto op:v.new_cigar){
    if(op.operation!='=')Emit(filtered,op.operation,op.length);
    else{
      std::vector<Range> excluded;
      for(auto [a,b]:uncertain){const auto lo=std::max(a,qp),hi=std::min(b,qp+op.length);if(lo<hi)excluded.emplace_back(lo,hi);}
      for(const auto&keep:*fixed){
        if(keep.q>=qp+op.length||keep.q+keep.n<=qp)continue;
        const auto lo=std::max(keep.q,qp),hi=std::min(keep.q+keep.n,qp+op.length);
        if(rp+lo-qp!=keep.r+lo-keep.q)continue;
        std::vector<Range> remaining;
        for(auto [a,b]:excluded){if(a<lo)remaining.emplace_back(a,std::min(b,lo));if(b>hi)remaining.emplace_back(std::max(a,hi),b);}
        excluded=std::move(remaining);
      }
      Position end=qp;
      for(auto [a,b]:Union(std::move(excluded))){Emit(filtered,'=',a-end);Emit(filtered,'D',b-a);Emit(filtered,'I',b-a);end=b;}
      Emit(filtered,'=',qp+op.length-end);
    }
    if(op.operation!='I')rp+=op.length;
    if(op.operation!='D')qp+=op.length;
  }
  std::vector<CigarOp> result;Length dr=0,dq=0;
  auto flush=[&]{Emit(result,'D',dr);Emit(result,'I',dq);dr=dq=0;};
  for(auto op:filtered){if(op.operation=='='){flush();Emit(result,'=',op.length);}else if(op.operation=='D')dr+=op.length;else if(op.operation=='I')dq+=op.length;else throw AlignmentError("nonexact retained component");}
  flush();v.new_cigar=std::move(result);
}

// Equal short opposing gaps provide a substitution explanation, not new
// evidence for locating an unknown repeated component. The zero-X route must
// not evade that ambiguity merely by spelling X as an equal I and D.
bool JointNewSubstitutionCorner(const CentromereEndpointReview&v,
    const CentromereStitchResult&old){
  Position rp=v.reference_begin,qp=v.query_begin;
  for(std::size_t i=0;i<v.new_cigar.size();++i){
    const auto a=v.new_cigar[i];
    if(i+1<v.new_cigar.size()){
      const auto b=v.new_cigar[i+1];
      if(((a.operation=='I'&&b.operation=='D')||(a.operation=='D'&&b.operation=='I'))&&a.length==b.length&&a.length<64){
        bool observed=false;
        for(const auto&segment:old.segments){
          const auto&path=segment.alignment;
          if(segment.path_role!="primary"||path.strand!=Strand::Forward)continue;
          Position rr=path.reference_begin,qq=path.query_begin;
          for(std::size_t j=0;j<path.cigar.size();++j){
            const auto op=path.cigar[j];
            if(op.operation=='X'&&rp>=rr&&qp>=qq&&rp-rr==qp-qq&&rp-rr+a.length<=op.length)observed=true;
            if(j+1<path.cigar.size()){
              const auto next=path.cigar[j+1];
              if(rr==rp&&qq==qp&&((op.operation=='I'&&next.operation=='D')||(op.operation=='D'&&next.operation=='I'))&&op.length==a.length&&next.length==a.length)observed=true;
            }
            if(op.operation!='I')rr+=op.length;
            if(op.operation!='D')qq+=op.length;
          }
        }
        if(!observed)return true;
      }
    }
    if(a.operation!='I')rp+=a.length;
    if(a.operation!='D')qp+=a.length;
  }
  return false;
}
// A maximal exact word can cross the boundary between overlapping physical
// copies. Ordered placement of the word does not assign those overlapping
// bases to an ancestral copy. Withhold only new overlap relations; retain old
// pairs and let the separate copy layer represent additional instances.
void JointMaskOverlappingCopies(CentromereEndpointReview&v,const SequenceRecord&r,
    const SequenceRecord&q,const CentromereStitchResult&old,
    const CentromereStitchOptions&options,bool withhold_over_budget){
  std::vector<Range> uncertain;Position rp=v.reference_begin,qp=v.query_begin;
  for(std::size_t oi=0;oi<v.new_cigar.size();++oi){
    const auto op=v.new_cigar[oi];
    if(op.operation=='='&&op.length>=64){
      const auto word=std::string_view(r.bases).substr(rp,op.length);
      auto overlaps=[&](const std::string&text,Position selected,bool query_side){
        auto at=text.find(word);unsigned count=0;
        while(at!=std::string::npos){
          if(options.interruption_callback)options.interruption_callback("joint-copy-overlap");
          if(++count>100){
            if(!withhold_over_budget)throw detail::StitchGapResourceLimit("context_occurrence_budget");
            uncertain.emplace_back(qp,qp+op.length);
            return count;
          }
          if(at!=selected){
            const auto lo=std::max<Position>(selected,at),hi=std::min<Position>(selected+op.length,at+op.length);
            if(lo<hi)uncertain.emplace_back(query_side?lo:qp+lo-rp,query_side?hi:qp+hi-rp);
          }
          at=text.find(word,at+1);
        }
        return count;
      };
      const auto nr=overlaps(r.bases,rp,false),nq=overlaps(q.bases,qp,true);
      auto opposing=[&](std::size_t i){
        if(i+1>=v.new_cigar.size())return false;
        const auto x=v.new_cigar[i],y=v.new_cigar[i+1];
        return ((x.operation=='I'&&y.operation=='D')||(x.operation=='D'&&y.operation=='I'))&&x.length==y.length;
      };
      if(nr<=100 && nq<=100 && nr!=nq && ((oi>=2&&opposing(oi-2))||opposing(oi+1))){
        auto counts=[&](std::string_view text,std::string_view needle){
          unsigned n=0;auto at=text.find(needle);
          while(at!=std::string::npos){
            if(options.interruption_callback)options.interruption_callback("joint-copy-context-core");
            if(++n>100){
              if(!withhold_over_budget)throw detail::StitchGapResourceLimit("context_occurrence_budget");
              // Both full-word counts are bounded, so this prefix cannot
              // belong to the same occurrence class.
              return n;
            }
            at=text.find(needle,at+1);
          }
          return n;
        };
        auto same=[&](std::string_view text){return counts(r.bases,text)==nr&&counts(q.bases,text)==nq;};
        auto minimum=[&](bool suffix){
          Length lo=32,hi=op.length;
          while(lo<hi){const auto mid=lo+(hi-lo)/2;
            if(same(suffix?word.substr(op.length-mid):word.substr(0,mid)))hi=mid;else lo=mid+1;
          }
          return lo;
        };
        const auto prefix=minimum(false),suffix=minimum(true);
        if(prefix>op.length-suffix){
          const auto first=op.length-suffix,last=prefix;
          // Opposing gaps plus unequal occurrence classes can create a
          // chimeric maximal extension across a copy boundary. Keep only the
          // occurrence-preserving core supported from both ends, or retain
          // ambiguity if that common core does not preserve the class.
          if(last-first>=32&&same(word.substr(first,last-first))){
            if(first)uncertain.emplace_back(qp,qp+first);
            if(last<op.length)uncertain.emplace_back(qp+last,qp+op.length);
          }else uncertain.emplace_back(qp,qp+op.length);
        }
      }
    }
    if(op.operation!='I')rp+=op.length;
    if(op.operation!='D')qp+=op.length;
  }
  if(uncertain.empty())return;
  uncertain=Union(std::move(uncertain));std::vector<CigarOp> filtered;rp=v.reference_begin;qp=v.query_begin;
  for(auto op:v.new_cigar){
    if(op.operation!='=')Emit(filtered,op.operation,op.length);
    else{
      std::vector<Range> excluded;
      for(auto [a,b]:uncertain){const auto lo=std::max(a,qp),hi=std::min(b,qp+op.length);if(lo<hi)excluded.emplace_back(lo,hi);}
      for(const auto&segment:old.segments){
        const auto&a=segment.alignment;if(segment.path_role!="primary"||a.strand!=Strand::Forward)continue;
        Position rr=a.reference_begin,qq=a.query_begin;
        for(auto previous:a.cigar){
          if(previous.operation=='='){
            const auto lo=std::max(qq,qp),hi=std::min(qq+previous.length,qp+op.length);
            if(lo<hi&&rp+lo-qp==rr+lo-qq){
              std::vector<Range> remaining;
              for(auto [x,y]:excluded){if(x<lo)remaining.emplace_back(x,std::min(y,lo));if(y>hi)remaining.emplace_back(std::max(x,hi),y);}
              excluded=std::move(remaining);
            }
          }
          if(previous.operation!='I')rr+=previous.length;
          if(previous.operation!='D')qq+=previous.length;
        }
      }
      Position end=qp;
      for(auto [a,b]:Union(std::move(excluded))){Emit(filtered,'=',a-end);Emit(filtered,'D',b-a);Emit(filtered,'I',b-a);end=b;}
      Emit(filtered,'=',qp+op.length-end);
    }
    if(op.operation!='I')rp+=op.length;
    if(op.operation!='D')qp+=op.length;
  }
  std::vector<CigarOp> result;Length dr=0,dq=0;
  auto flush=[&]{Emit(result,'D',dr);Emit(result,'I',dq);dr=dq=0;};
  for(auto op:filtered){if(op.operation=='='){flush();Emit(result,'=',op.length);}else if(op.operation=='D')dr+=op.length;else if(op.operation=='I')dq+=op.length;else throw AlignmentError("nonexact copy overlap path");}
  flush();v.new_cigar=std::move(result);
}

// A repeat word need not locate itself. Require a unique ordered assignment of
// multiple independently interrupted words between the same physical flanks.
// Every occurrence is retained within the bounded problem; exhausting a limit
// is unresolved, never evidence of uniqueness. This certifies these contexts,
// not unrestricted global alignment optimality.
bool JointOrderedContexts(const SequenceRecord&r,const SequenceRecord&q,
    const CentromereEndpointReview&v,const CentromereStitchOptions&o){
  struct State {Position r{},q{};std::int64_t cost{};unsigned ways{};bool selected{};};
  std::vector<State> prior{{v.reference_begin,v.query_begin,0,1,true}};
  Length previous_length=0;Position rp=v.reference_begin,qp=v.query_begin;
  unsigned layers=0;std::size_t visits=0,nodes=0;std::set<std::string> words;
  auto penalty=[](Length n)->std::int64_t{return n?std::min({60+30*static_cast<std::int64_t>(n),800+5*static_cast<std::int64_t>(n),2500+static_cast<std::int64_t>(n)}):0;};
  auto gap=[&](Length a,Length b){return penalty(a>b?a-b:b-a);};
  auto occurrences=[&](const std::string&text,const std::string&word,Position lo,Position hi){
    std::vector<Position> hits;unsigned count=0;std::size_t at=text.find(word);
    while(at!=std::string::npos){
      if(o.interruption_callback)o.interruption_callback("joint-ordered-occurrences");
      if(++count>100)throw detail::StitchGapResourceLimit("context_occurrence_budget");
      if(at>=lo&&at+word.size()<=hi)hits.push_back(at);
      at=text.find(word,at+1);
    }
    return hits;
  };
  for(auto op:v.new_cigar){
    if(op.operation=='='&&op.length>=64){
      auto word=r.bases.substr(rp,op.length);
      // Identical or nested words are not independent witnesses.
      bool redundant=false;for(const auto&old:words)if(old.find(word)!=std::string::npos||word.find(old)!=std::string::npos)redundant=true;
      if(!redundant){
        auto rh=occurrences(r.bases,word,v.reference_begin,v.reference_end);
        auto qh=occurrences(q.bases,word,v.query_begin,v.query_end);
        std::vector<State> next;
        for(auto rr:rh)for(auto qq:qh){
          if(++nodes>3000)throw detail::StitchGapResourceLimit("candidate_budget");
          State x{rr,qq,std::numeric_limits<std::int64_t>::max(),0,false};
          for(const auto&old:prior){
            if(++visits>100000)throw detail::StitchGapResourceLimit("candidate_budget");
            if(old.r+previous_length>rr||old.q+previous_length>qq||!old.ways)continue;
            const auto cost=old.cost+gap(rr-old.r-previous_length,qq-old.q-previous_length);
            const bool selected=old.selected&&rr==rp&&qq==qp;
            if(cost<x.cost){x.cost=cost;x.ways=old.ways;x.selected=selected;}
            else if(cost==x.cost){x.ways=std::min(2U,x.ways+old.ways);x.selected=x.selected||selected;}
          }
          if(x.ways)next.push_back(x);
        }
        if(next.empty())return false;
        prior=std::move(next);previous_length=op.length;words.insert(std::move(word));++layers;
      }
    }
    if(op.operation!='I')rp+=op.length;
    if(op.operation!='D')qp+=op.length;
  }
  if(layers<2)return false;
  std::int64_t best=std::numeric_limits<std::int64_t>::max();unsigned ways=0;bool selected=false;
  for(const auto&x:prior){
    const auto cost=x.cost+gap(v.reference_end-x.r-previous_length,v.query_end-x.q-previous_length);
    if(cost<best){best=cost;ways=x.ways;selected=x.selected;}
    else if(cost==best){ways=std::min(2U,ways+x.ways);selected=selected||x.selected;}
  }
  return ways==1&&selected;
}

// A global chain can be ordered while a particular repeat word is on the
// wrong copy. Certify its two locally identifying ends independently, using
// the occurrence ranks outside the proposed corridor as a consistency check.
bool JointUnplannedContext(const SequenceRecord&r,const SequenceRecord&q,
    const CentromereEndpointReview&v,JointSeed seed,const CentromereStitchOptions&options){
  const auto reference=std::string_view(r.bases);
  const auto query=std::string_view(q.bases);
  const auto rview=reference.substr(v.reference_begin,v.reference_end-v.reference_begin);
  const auto qview=query.substr(v.query_begin,v.query_end-v.query_begin);
  const auto word=reference.substr(seed.r,seed.n);
  auto unique=[](std::string_view text,std::string_view pattern){
    const auto at=text.find(pattern);
    return at!=std::string_view::npos && text.find(pattern,at+1)==std::string_view::npos;
  };
  if(!unique(rview,word) || !unique(qview,word))return false;
  auto shortest=[&](bool suffix){
    Length lo=32,hi=seed.n;
    while(lo<hi){
      if(options.interruption_callback)options.interruption_callback("joint-unplanned-context");
      const auto mid=lo+(hi-lo)/2;
      const auto pattern=suffix?word.substr(seed.n-mid):word.substr(0,mid);
      if(unique(rview,pattern) && unique(qview,pattern))hi=mid;
      else lo=mid+1;
    }
    return lo;
  };
  const auto prefix=shortest(false),suffix=shortest(true);
  // The locators must not share any full 32-base evidence word.
  if(prefix+suffix>=seed.n+32)return false;
  auto signature=[&](std::string_view text,std::string_view pattern,Position selected){
    unsigned count=0,rank=0;
    for(auto at=text.find(pattern);at!=std::string_view::npos;at=text.find(pattern,at+1)){
      if(options.interruption_callback)options.interruption_callback("joint-unplanned-context-rank");
      if(++count>100)break;
      rank+=at<selected;
    }
    return std::pair{count,rank};
  };
  auto consistent=[&](Length offset,Length length){
    const auto pattern=word.substr(offset,length);
    const auto a=signature(reference,pattern,seed.r+offset);
    return a.first>0 && a.first<=100 && a==signature(query,pattern,seed.q+offset);
  };
  return consistent(0,prefix) && consistent(seed.n-suffix,suffix);
}

// Project only exact, unclaimed coordinates from an ordered path into a
// direction-unknown interval. The surrounding path is a localization witness;
// it does not authorize replacing the existing segments in that rectangle.
bool JointSelectUnplannedRun(CentromereEndpointReview&v,
    const CentromereDirectedAnchorResult&anchors,const CentromereStitchResult&old,
    const SequenceRecord&r,const SequenceRecord&q,const CentromereStitchOptions&options){
  JointSeed best{};
  Position rp=v.reference_begin,qp=v.query_begin;
  for(auto op:v.new_cigar){
    if(op.operation=='=' && op.length>=64 &&
       JointUnplannedContext(r,q,v,{rp,qp,op.length,0,0},options))
      for(const auto&unknown:anchors.unknown){
      const auto qb=std::max(qp,unknown.query_begin),qe=std::min(qp+op.length,unknown.query_end);
      if(qb>=qe || qe-qb<64 || qe-qb<best.n)continue;
      std::vector<Range> occupied;
      const auto rb=rp+qb-qp,re=rb+qe-qb;
      for(const auto&segment:old.segments){
        if(options.interruption_callback)options.interruption_callback("joint-unplanned-footprint");
        const auto&a=segment.alignment;
        // Reserve actual old pair coordinates. A long deletion can span this
        // reference interval while pairing none of its bases.
        Position rr=a.reference_begin;
        for(const auto&old_op:a.cigar){
          if((old_op.operation=='=' || old_op.operation=='X') && Touch(rb,re,rr,rr+old_op.length))
            occupied.emplace_back(qb+std::max(rb,rr)-rb,qb+std::min(re,rr+old_op.length)-rb);
          if(old_op.operation!='I')rr+=old_op.length;
        }
        if(Touch(qb,qe,a.query_begin,a.query_end))
          occupied.emplace_back(std::max(qb,a.query_begin),std::min(qe,a.query_end));
      }
      auto add=[&](Position begin,Position end){
        const auto n=end-begin;
        if(n<64 || n<=best.n)return;
        const auto rr=rp+begin-qp;
        const auto word=std::string_view(r.bases).substr(rr,n);
        auto occurrences=[&](const std::string&sequence){
          unsigned count=0;
          for(auto at=sequence.find(word);at!=std::string::npos;at=sequence.find(word,at+1)){
            if(options.interruption_callback)options.interruption_callback("joint-unplanned-occurrences");
            if(++count>100)break;
          }
          return count;
        };
        const auto nr=occurrences(r.bases);
        if(nr>0 && nr<=100 && occurrences(q.bases)==nr)best={rr,begin,n,nr,nr};
      };
      auto begin=qb;
      for(const auto&[lo,hi]:Union(std::move(occupied))){add(begin,lo);begin=hi;}
      add(begin,qe);
    }
    if(op.operation!='I')rp+=op.length;
    if(op.operation!='D')qp+=op.length;
  }
  if(!best.n)return false;
  v.reference_begin=best.r;v.reference_end=best.r+best.n;
  v.query_begin=best.q;v.query_end=best.q+best.n;
  v.new_cigar={{'=',best.n}};
  return true;
}

CentromereDirectedAnchorResult JointReplace(const CentromereDirectedAnchorResult &old,
    const CentromereEndpointReview &v,const SequenceRecord &r,const SequenceRecord &q){
  auto result=old;
  auto allocate=[](std::uint64_t&value){
    if(value==std::numeric_limits<std::uint64_t>::max())throw AlignmentError("joint review identity overflow");
    return ++value;
  };
  std::uint64_t next=0;
  for(const auto &region:old.regions)next=std::max(next,region.id);
  const auto added=allocate(next);
  result.regions.clear();
  for(auto region:old.regions){
    if(!Touch(region.query_begin,region.query_end,v.query_begin,v.query_end)){
      result.regions.push_back(region);continue;
    }
    const auto finish=region.query_end;
    if(region.query_begin<v.query_begin){region.query_end=v.query_begin;result.regions.push_back(region);}
    if(finish>v.query_end){region.id=allocate(next);region.query_begin=v.query_end;region.query_end=finish;result.regions.push_back(region);}
  }
  result.regions.push_back({added,v.query_begin,v.query_end,Strand::Forward,0});
  std::sort(result.regions.begin(),result.regions.end(),[](const auto&a,const auto&b){return a.query_begin<b.query_begin;});
  auto replace=[&](std::vector<A> &list){
    std::uint64_t id=0;for(const auto &a:list)id=std::max(id,a.id);
    std::vector<A> output;
    for(auto a:list){
      if(Touch(a.query_begin,a.query_begin+a.length,v.query_begin,v.query_end))continue;
      bool found=false;
      for(const auto &region:result.regions)
        if(a.query_begin>=region.query_begin && a.query_begin+a.length<=region.query_end){a.region_id=region.id;found=true;break;}
      if(!found)throw AlignmentError("joint exterior anchor lost its direction region");
      output.push_back(a);
    }
    Position rp=v.reference_begin,qp=v.query_begin;
    for(const auto op:v.new_cigar){
      if(op.operation=='=')output.push_back({allocate(id),added,0,rp,qp,op.length,Strand::Forward,static_cast<double>(op.length)});
      if(op.operation!='I')rp+=op.length;
      if(op.operation!='D')qp+=op.length;
    }
    std::sort(output.begin(),output.end(),[](const A&a,const A&b){return std::tuple{a.region_id,a.reference_begin,a.query_begin}<std::tuple{b.region_id,b.reference_begin,b.query_begin};});
    std::map<std::uint64_t,std::uint64_t> rank;
    for(auto &a:output){a.chain_index=++rank[a.region_id];CheckAnchor(a,r,q);}
    list=std::move(output);
  };
  replace(result.initial);replace(result.final);
  result.tasks.clear();
  std::set<std::uint64_t> retained;
  for(const auto &region:result.regions){
    retained.insert(region.id);
    CentromereDirectedTask task;
    for(const auto &prior:old.tasks)if(prior.region.id==region.id){task=prior;break;}
    task.region=region;task.reference_length=r.size();
    task.reference_begin=region.id==added?v.reference_begin:0;
    task.reference_end=region.id==added?v.reference_end:r.size();
    task.status="complete";result.tasks.push_back(task);
  }
  std::erase_if(result.trace,[&](const auto &x){
    if(!retained.contains(x.region_id))return true;
    const auto prior=std::find_if(old.regions.begin(),old.regions.end(),[&](const auto&a){return a.id==x.region_id;});
    const auto current=std::find_if(result.regions.begin(),result.regions.end(),[&](const auto&a){return a.id==x.region_id;});
    return prior==old.regions.end() || current==result.regions.end() || *prior!=*current;
  });
  for(auto&task:result.tasks){
    task.recursion_node_count=0;task.max_recursion_depth=0;
    for(const auto&trace:result.trace)if(trace.region_id==task.region.id){++task.recursion_node_count;task.max_recursion_depth=std::max(task.max_recursion_depth,trace.node.depth);}
  }
  result.unknown.clear();Position end=0;
  for(const auto &region:result.regions){if(end<region.query_begin)result.unknown.push_back({0,end,region.query_begin,Strand::Forward});end=region.query_end;}
  if(end<q.size())result.unknown.push_back({0,end,q.size(),Strand::Forward});
  Refresh(result);return result;
}
std::vector<CentromereStitchGap> JointReplay(const CentromereStitchResult &old,
    const CentromereDirectedAnchorResult &anchors,const SequenceRecord&r,
    const SequenceRecord&q,const CentromereEndpointReview&v){
  auto result=CrossFixedGaps(r,q,anchors,v);
  std::map<std::uint64_t,std::vector<const A*>> paths;
  for(const auto&a:anchors.final)paths[a.region_id].push_back(&a);
  for(auto &[region,list]:paths){
    std::sort(list.begin(),list.end(),[](const A*a,const A*b){return a->chain_index<b->chain_index;});
    for(std::size_t i=1;i<list.size();++i){
      const auto&a=*list[i-1],&b=*list[i];
      const Position rb=a.reference_begin+a.length,re=b.reference_begin;
      const Position qb=a.strand==Strand::Forward?a.query_begin+a.length:b.query_begin+b.length;
      const Position qe=a.strand==Strand::Forward?b.query_begin:a.query_begin;
      if(rb==re&&qb==qe)continue;
      if(rb>=v.reference_begin&&re<=v.reference_end&&qb>=v.query_begin&&qe<=v.query_end)continue;
      for(auto gap:old.gaps)
        if(gap.path_role=="primary"&&gap.strand==a.strand&&gap.reference_begin==rb&&gap.reference_end==re&&gap.query_begin==qb&&gap.query_end==qe&&gap.status=="complete"&&gap.scored&&
           std::none_of(gap.cigar.begin(),gap.cigar.end(),[](auto op){return op.operation=='X';})){
          gap.region_id=region;result.push_back(std::move(gap));break;
        }
    }
  }
  for(const auto&gap:old.gaps)if(gap.path_role=="copy"&&gap.status=="complete"&&gap.scored&&
      std::none_of(gap.cigar.begin(),gap.cigar.end(),[](auto op){return op.operation=='X';}))result.push_back(gap);
  return result;
}
#endif
} // namespace

CentromereEndpointResult RefineCentromereJointPaths(
    const SequenceRecord&r,const SequenceRecord&q,const CentromereDirectedAnchorResult&input,
    const CentromereCopyAnchorResult*copies,const CentromereStitchResult&stitched,
    const CentromereEndpointOptions&options){
#if !RAMAG_HAVE_CENTROMERE
  (void)r;(void)q;(void)input;(void)copies;(void)stitched;(void)options;
  throw DependencyError("joint path review requires the source RaMA backend");
#else
  if(!input.complete || (copies&&!copies->complete) || stitched.aborted || !stitched.complete ||
     stitched.reference_length!=r.size() || stitched.query_length!=q.size())
    throw AlignmentError("joint path review requires complete valid stages");
  if(!options.max_reviews || !std::isfinite(options.stage_seconds_limit) || options.stage_seconds_limit<=0 ||
     !std::isfinite(options.window_seconds_limit)||options.window_seconds_limit<=0 ||
     !options.kernel.gap_memory_limit_bytes || !std::isfinite(options.kernel.gap_seconds_limit)||options.kernel.gap_seconds_limit<=0 ||
     options.kernel.score_model!=CentromereStitchScoreModel::Centrolign3Piece || stitched.score_model!=options.kernel.score_model)
    throw AlignmentError("invalid joint path scoring/resources");
  const auto started=Clock::now();CentromereEndpointResult out;out.anchors=input;
  auto check=[&]{if(options.kernel.interruption_callback)options.kernel.interruption_callback("centromere-joint-path-review");};
  const auto rc=CountContexts(r,options),qc=CountContexts(q,options);
  SequenceRecord reverse=q;reverse.bases=CrossReverseComplement(q.bases);
  const auto qrc=CountContexts(reverse,options);
  std::vector<A> flanks;
  for(const auto&a:input.final){CheckAnchor(a,r,q);if(a.strand==Strand::Forward && (DistinguishableContext(a,r,rc,qc)||JointUnique(a,r,q)))flanks.push_back(a);}
  std::sort(flanks.begin(),flanks.end(),[](const A&a,const A&b){return std::tuple{a.query_begin,a.reference_begin}<std::tuple{b.query_begin,b.reference_begin};});
  std::vector<Range> used_r,used_q;
  struct UnresolvedContext {Position rb{},re{},qb{},qe{};};
  std::vector<UnresolvedContext> unresolved_contexts, marker_conflicts, unplanned_contexts;
  std::map<std::pair<Position,Length>,std::string> multiplicity;
  auto copy_competition=[&](const A&a){
    if(a.strand!=Strand::Forward||a.length<64)return std::string{};
    const auto key=std::pair{a.reference_begin,a.length};
    if(const auto found=multiplicity.find(key);found!=multiplicity.end())return found->second;
    const auto word=std::string_view(r.bases).substr(a.reference_begin,a.length);
    auto count=[&](const std::string&text){
      std::uint64_t n=0;std::size_t start=0;
      while(n<=100){check();const auto at=text.find(word,start);if(at==std::string::npos)break;++n;start=at+1;}
      return n;
    };
    const auto nr=count(r.bases),nq=count(q.bases);
    if(!nr||!nq)throw AlignmentError("joint original multiplicity is zero");
    std::string reason;
    if(nr>100||nq>100)reason="primary_multiplicity_count_budget";
    else if(nr!=nq)reason="primary_multiplicity_requires_copy_review";
    multiplicity.emplace(key,reason);return reason;
  };
  struct Window {std::size_t left{},right{};Length unknown{};unsigned strong{};};
  std::vector<Window> windows;
  for(std::size_t hop=1;hop<=4;++hop)for(std::size_t i=0;i+hop<flanks.size();++i){
    const auto &a=flanks[i],&b=flanks[i+hop];
    if(a.query_begin+a.length>=b.query_begin || a.reference_begin+a.length>=b.reference_begin)continue;
    Length unknown=0;
    for(const auto &x:input.unknown){const auto lo=std::max(a.query_begin+a.length,x.query_begin),hi=std::min(b.query_begin,x.query_end);if(lo<hi)unknown+=hi-lo;}
    windows.push_back({i,i+hop,unknown,static_cast<unsigned>(DistinguishableContext(a,r,rc,qc))+static_cast<unsigned>(DistinguishableContext(b,r,rc,qc))});
  }
  // Weak maximal matches must not hide the nearest independently positioned
  // flanks simply by occupying all four bounded neighbor slots.
  std::vector<std::size_t> positioned;
  std::set<std::pair<std::size_t,std::size_t>> seen;
  for(const auto&v:windows)seen.emplace(v.left,v.right);
  for(std::size_t i=0;i<flanks.size();++i)if(DistinguishableContext(flanks[i],r,rc,qc))positioned.push_back(i);
  for(std::size_t hop=1;hop<=4;++hop)for(std::size_t j=0;j+hop<positioned.size();++j){
    const auto a=positioned[j],b=positioned[j+hop];
    if(!seen.emplace(a,b).second)continue;
    const auto &left=flanks[a],&right=flanks[b];
    if(left.query_begin+left.length>=right.query_begin || left.reference_begin+left.length>=right.reference_begin)continue;
    Length unknown=0;
    for(const auto&x:input.unknown){const auto lo=std::max(left.query_begin+left.length,x.query_begin),hi=std::min(right.query_begin,x.query_end);if(lo<hi)unknown+=hi-lo;}
    windows.push_back({a,b,unknown,static_cast<unsigned>(DistinguishableContext(left,r,rc,qc))+static_cast<unsigned>(DistinguishableContext(right,r,rc,qc))});
  }
  // Resolve a complete conflicting interval before its weak inner roots can
  // become fixed outer boundaries for smaller, overlapping reviews.
  std::stable_sort(windows.begin(),windows.end(),[&](const Window&a,const Window&b){
    if(a.unknown!=b.unknown)return a.unknown>b.unknown;
    if(a.strong!=b.strong)return a.strong>b.strong;
    const auto span=[&](const Window&v){return flanks[v.right].query_begin-flanks[v.left].query_begin-flanks[v.left].length;};
    return span(a)>span(b);
  });
  // Preserve the existing pass before trying stronger composite evidence.
  // Both passes share stage time, DP calls, footprints and accepted output.
  bool constrained=false;
  bool core_review=false;
  bool unplanned_review=false;
  std::uint64_t unplanned_attempts=0;
  for(std::uint64_t composite=0;composite<=options.max_reviews;++composite){
    bool changed=false;
    for(const auto &window:windows){
    check();
    if(Seconds(started)>=options.stage_seconds_limit ||
       (unplanned_review ? unplanned_attempts : out.kernel_calls)>=options.max_reviews)break;
    const auto&left=flanks[window.left];const auto&right=flanks[window.right];
    const Position rb=left.reference_begin+left.length,re=right.reference_begin,qb=left.query_begin+left.length,qe=right.query_begin;
    if(rb>=re||qb>=qe||re-rb>65536||qe-qb>65536)continue;
    if(unplanned_review && std::none_of(unplanned_contexts.begin(),unplanned_contexts.end(),[&](const auto&x){
      return rb==x.rb && re==x.re && qb==x.qb && qe==x.qe;
    }))continue;
    if(core_review && std::none_of(marker_conflicts.begin(),marker_conflicts.end(),[&](const auto&x){
      return rb==x.rb && re==x.re && qb==x.qb && qe==x.qe;
    }))continue;
    bool overlap=false;
    for(std::size_t j=0;j<used_r.size();++j)overlap|=Touch(rb,re,used_r[j].first,used_r[j].second)||Touch(qb,qe,used_q[j].first,used_q[j].second);
    if(overlap)continue;
    CentromereEndpointReview v;v.id=out.reviews.size()+1;v.region_id=left.region_id;v.left_anchor_id=left.id;v.right_anchor_id=right.id;v.reference_begin=rb;v.reference_end=re;v.query_begin=qb;v.query_end=qe;
    const auto begin=Clock::now();
    // A smaller rectangle cannot discard the discordant part of an enclosing
    // independently bounded context that contradicted a protected marker or failed relative consistency.
    // Keep that whole local competition unresolved in this pass.
    if(!constrained && !unplanned_review && std::any_of(unresolved_contexts.begin(),unresolved_contexts.end(),[&](const auto&x){
         return rb>=x.rb && re<=x.re && qb>=x.qb && qe<=x.qe;
       })){
      v.reason="enclosing_context_consistency_unresolved";out.reviews.push_back(std::move(v));continue;
    }
    const bool exterior=DistinguishableContext(left,r,rc,qc)&&DistinguishableContext(right,r,rc,qc);
    const Length displacement=(re-rb)>(qe-qb)?(re-rb)-(qe-qb):(qe-qb)-(re-rb);
    if(!exterior && displacement>=64 && !composite){
      v.reason="no_independent_exterior_context";out.reviews.push_back(std::move(v));continue;
    }
    if(composite && (!JointVariableFlank(left,r,q)||!JointVariableFlank(right,r,q)))continue;
    const auto &current=out.anchors;const auto &assembly=out.stitched?*out.stitched:stitched;
    bool conflict=false;std::size_t inside=0;
    for(const auto&a:current.final){
      const bool tq=Touch(a.query_begin,a.query_begin+a.length,qb,qe),tr=Touch(a.reference_begin,a.reference_begin+a.length,rb,re);
      if(!tq&&!tr)continue;
      if(!tq || a.query_begin<qb || a.query_begin+a.length>qe){
        if(!unplanned_review){v.reason="exterior_anchor_crosses_window";break;}
      }
      ++inside;
      const bool strong=a.strand==Strand::Forward?DistinguishableContext(a,r,rc,qc):DistinguishableContext(a,r,rc,qrc);
      if(a.strand==Strand::Reverse){conflict=true;auto oriented=a;oriented.query_begin=q.size()-a.query_begin-a.length;
        if(strong||JointUnique(oriented,r,reverse)){v.reason="protected_reverse_context";break;}}
      if((a.reference_begin<rb||a.reference_begin+a.length>re)&&strong){v.reason="protected_exterior_context";break;}
    }
    const bool has_unknown=std::any_of(current.unknown.begin(),current.unknown.end(),
      [&](const auto&x){return Touch(qb,qe,x.query_begin,x.query_end);});
    if(unplanned_review && !has_unknown)continue;
    conflict|=has_unknown;
    if(has_unknown && v.reason=="exterior_anchor_crosses_window")
      unplanned_contexts.push_back({rb,re,qb,qe});
    std::vector<const A*> forward_context;
    for(const auto&a:current.final)if(a.strand==Strand::Forward&&a.query_begin>=qb&&a.query_begin+a.length<=qe)forward_context.push_back(&a);
    std::sort(forward_context.begin(),forward_context.end(),[](const A*a,const A*b){return a->query_begin<b->query_begin;});
    for(std::size_t j=1;j<forward_context.size();++j)if(forward_context[j-1]->reference_begin+forward_context[j-1]->length>forward_context[j]->reference_begin)conflict=true;
    if(inside>options.kernel.gap_memory_limit_bytes/(8*sizeof(A)))v.reason="block_geometry_budget";
    if(!conflict)continue;
    std::vector<A> protected_copy_roots;
    if(v.reason.empty())for(const auto&a:current.initial)
      if(Touch(a.query_begin,a.query_begin+a.length,qb,qe)||Touch(a.reference_begin,a.reference_begin+a.length,rb,re)){
        const auto reason=copy_competition(a);
        if(reason=="primary_multiplicity_requires_copy_review")protected_copy_roots.push_back(a);
        else if(!reason.empty()){v.reason=reason;break;}
      }
    if(copies)for(const auto&a:copies->anchors)if(Touch(rb,re,a.reference_begin,a.reference_begin+a.length)||Touch(qb,qe,a.query_begin,a.query_begin+a.length))v.reason="accepted_copy_footprint";
    for(const auto&s:assembly.segments)if(s.path_role=="copy"&&(Touch(rb,re,s.alignment.reference_begin,s.alignment.reference_end)||Touch(qb,qe,s.alignment.query_begin,s.alignment.query_end)))v.reason="accepted_copy_segment_footprint";
    if(!v.reason.empty()){out.reviews.push_back(std::move(v));continue;}
    auto old=CrossFragments(assembly,r,q,rb,re,qb,qe);
    auto reversed=assembly;
    std::erase_if(reversed.segments,[](const auto&s){return s.path_role!="primary"||s.alignment.strand!=Strand::Reverse;});
    for(auto&s:reversed.segments){auto&a=s.alignment;const auto b=a.query_begin;a.query_begin=q.size()-a.query_end;a.query_end=q.size()-b;a.strand=Strand::Forward;}
    const auto old_reverse=CrossFragments(reversed,r,reverse,rb,re,q.size()-qe,q.size()-qb);
    old.score+=old_reverse.score;old.exact+=old_reverse.exact;old.mismatches+=old_reverse.mismatches;old.short_disagreements+=old_reverse.short_disagreements;
    v.old_score=old.score;v.old_exact_bases=old.exact;
    if(old.exact>=std::min(re-rb,qe-qb))continue;
    if(r.bases.compare(rb,re-rb,q.bases,qb,qe-qb)==0)continue;
    auto kernel=options.kernel;kernel.observation_callback={};kernel.progress_callback={};
    const auto budget=std::min(options.window_seconds_limit,options.stage_seconds_limit-Seconds(started));
    kernel.gap_seconds_limit=std::min(kernel.gap_seconds_limit,budget);kernel.stage_seconds_limit=budget;
    const auto caller=kernel.interruption_callback;
    kernel.interruption_callback=[&](std::string_view phase){if(caller)caller(phase);if(Seconds(begin)>budget)throw detail::StitchGapResourceLimit("gap_time_budget");};
    try{
      if(unplanned_review)++unplanned_attempts;
      const auto prior_contexts=CrossOldContexts(r,rc,qc,assembly,v);
      auto reciprocal=[&](Position part_rb,Position part_re,Position part_qb,Position part_qe){
        std::vector<JointSeed> consensus;
        if(part_rb==part_re || part_qb==part_qe)return consensus;
        CentromereEndpointReview part;
        part.reference_begin=part_rb;part.reference_end=part_re;part.query_begin=part_qb;part.query_end=part_qe;
        const auto seeds=JointSeeds(r,q,part,kernel);
        if(seeds.empty())return consensus;
        SequenceRecord mirrored_r,mirrored_q;
        mirrored_r.bases=CrossReverseComplement(r.bases.substr(part_rb,part_re-part_rb));
        mirrored_q.bases=CrossReverseComplement(q.bases.substr(part_qb,part_qe-part_qb));
        CentromereEndpointReview mirrored;mirrored.reference_end=part_re-part_rb;mirrored.query_end=part_qe-part_qb;
        auto counter=JointSeeds(mirrored_r,mirrored_q,mirrored,kernel);
        // Keep only seed correspondences endorsed from both traversal ends.
        // Distinct clipping choices are ambiguity, not a reason to force either
        // maximal-match endpoint or discard the whole informative window.
        for(auto &b:counter){b.r=part_re-b.r-b.n;b.q=part_qe-b.q-b.n;
          for(const auto&a:seeds){
            const auto lo=std::max(a.q,b.q),hi=std::min(a.q+a.n,b.q+b.n);
            if(lo<hi && a.r+lo-a.q==b.r+lo-b.q && hi-lo>=32)
              consensus.push_back({a.r+lo-a.q,lo,hi-lo,a.rc,a.qc});
          }
        }
        std::sort(consensus.begin(),consensus.end(),[](const auto&a,const auto&b){return std::tuple{a.r,a.q}<std::tuple{b.r,b.q};});
        return consensus;
      };
      std::vector<JointSeed> consensus;
      if(core_review || unplanned_review){
        // Keep actual independently locating pairs at their old coordinates.
        // Merge overlapping words into one physical core; count them once.
        std::vector<JointSeed> cores;
        for(const auto&[rr,qq]:prior_contexts)cores.push_back({rr,qq,32,1,1});
        const bool collinear=JointMergeCollinear(cores);
        if(!collinear || (core_review && cores.empty())){
          v.reason="old_cores_not_collinear";out.reviews.push_back(std::move(v));continue;
        }
        Position rr=rb,qq=qb;
        for(const auto&core:cores){
          const auto seeds=reciprocal(rr,core.r,qq,core.q);
          consensus.insert(consensus.end(),seeds.begin(),seeds.end());
          consensus.push_back(core);
          rr=core.r+core.n;qq=core.q+core.n;
        }
        const auto seeds=reciprocal(rr,re,qq,qe);
        consensus.insert(consensus.end(),seeds.begin(),seeds.end());
      }else consensus=reciprocal(rb,re,qb,qe);
      if(consensus.empty()){v.reason="no_reciprocal_rare_context";v.seconds=Seconds(begin);out.reviews.push_back(std::move(v));continue;}
      if(constrained){
        const auto fixed=JointConstrainPairs(assembly,v,consensus);
        if(!fixed){v.reason="old_pair_constraints_not_collinear";out.reviews.push_back(std::move(v));continue;}
        consensus=*fixed;
      }
      auto located_vertex=[&](Position rr,Position qq,bool left_side){
        if((rr==rb&&qq==qb)||(rr==re&&qq==qe))return true;
        for(const auto&a:consensus){
          if(left_side?(a.r+a.n!=rr||a.q+a.n!=qq):(a.r!=rr||a.q!=qq))continue;
          if(a.rc&&a.qc)return true;
          A anchor;anchor.reference_begin=a.r;anchor.query_begin=a.q;anchor.length=a.n;anchor.strand=Strand::Forward;
          if(DistinguishableContext(anchor,r,rc,qc)||JointVariableFlank(anchor,r,q))return true;
        }
        return false;
      };
      Position rp=rb,qp=qb;Length supported=0;
      std::vector<Range> confirmed;
      auto solve=[&](Position rend,Position qend){
        if(rend<rp||qend<qp)throw AlignmentError("joint seed order");
        // Unplanned coverage uses only exact seed pairs supported from both
        // ends. Do not run DP in the unevidenced space between them.
        if(unplanned_review){Emit(v.new_cigar,'D',rend-rp);Emit(v.new_cigar,'I',qend-qp);}
        else if(rend==rp)Emit(v.new_cigar,'I',qend-qp);
        else if(qend==qp)Emit(v.new_cigar,'D',rend-rp);
        else{
          // An unchanged, already verified interval is not a new search. Its
          // actual coordinate vertices and operation order must match exactly.
          // Replaying it preserves old repeat choices and reserves the shared
          // budget for the conflicting part of a larger context window.
          for(const auto&segment:assembly.segments){
            if(segment.path_role!="primary"||segment.alignment.strand!=Strand::Forward)continue;
            std::vector<CigarOp> previous;
            if(!SliceStitchedCigar(assembly,segment.region_id,rp,rend,qp,qend,previous)||
               std::any_of(previous.begin(),previous.end(),[](auto op){return op.operation=='X';}))continue;
            (void)detail::ScoreCentromereStitchCigar(std::string_view(r.bases).substr(rp,rend-rp),std::string_view(q.bases).substr(qp,qend-qp),previous,kernel.score_model);
            Position pr=rp,pq=qp;
            for(auto op:previous){Emit(v.new_cigar,op.operation,op.length);if(op.operation=='=')confirmed.emplace_back(pq,pq+op.length);if(op.operation!='I')pr+=op.length;if(op.operation!='D')pq+=op.length;}
            if(pr!=rend||pq!=qend)throw AlignmentError("joint replay consumption differs");
            rp=rend;qp=qend;return;
          }
          // An old short exact hit may be retained without locating a new
          // extension. New gap pairs require independently supported vertices
          // on both sides; unresolved tails keep only their consumption.
          if(constrained && (!located_vertex(rp,qp,true)||!located_vertex(rend,qend,false))){
            Emit(v.new_cigar,'D',rend-rp);Emit(v.new_cigar,'I',qend-qp);
            rp=rend;qp=qend;return;
          }
          if(out.kernel_calls>=options.max_reviews)throw detail::StitchGapResourceLimit("review_budget");
          ++out.kernel_calls;
          auto gap=detail::AlignCentromereBlockDpExactGaps(std::string_view(r.bases).substr(rp,rend-rp),std::string_view(q.bases).substr(qp,qend-qp),kernel);
          v.peak_kernel_bytes=std::max(v.peak_kernel_bytes,gap.peak_kernel_bytes);v.dp_cells+=gap.dp_cells;
          const bool opposing=[&]{for(std::size_t k=1;k<gap.cigar.size();++k)if((gap.cigar[k-1].operation=='I'&&gap.cigar[k].operation=='D')||(gap.cigar[k-1].operation=='D'&&gap.cigar[k].operation=='I'))return true;return false;}();
          // At independently located vertices the exact-gap and unrestricted
          // traces below must agree on new pairs. The short counterevidence
          // remains unpaired; it cannot establish a new locus.
          if(core_review && opposing &&
             (!located_vertex(rp,qp,true) || !located_vertex(rend,qend,false))){
            CentromereEndpointReview component;
            component.reference_begin=rp;component.reference_end=rend;
            component.query_begin=qp;component.query_end=qend;
            component.new_cigar=gap.cigar;
            if(JointNewSubstitutionCorner(component,assembly)){
              // This component has counterevidence. Keep its old physical
              // pairs, if they are collinear, without adopting new ones.
              const auto prior=CrossFragments(assembly,r,q,rp,rend,qp,qend);
              const auto fixed=JointConstrainPairs(assembly,component,{});
              if(prior.mismatches || (prior.exact && !fixed))
                throw detail::StitchGapResourceLimit("internal_core_component_unresolved");
              if(fixed)for(const auto&seed:*fixed){
                Emit(v.new_cigar,'D',seed.r-rp);Emit(v.new_cigar,'I',seed.q-qp);
                Emit(v.new_cigar,'=',seed.n);
                rp=seed.r+seed.n;qp=seed.q+seed.n;
              }
              Emit(v.new_cigar,'D',rend-rp);Emit(v.new_cigar,'I',qend-qp);
              rp=rend;qp=qend;return;
            }
          }
          std::vector<JointSeed> peer;
          if(opposing){
            if(out.kernel_calls>=options.max_reviews)throw detail::StitchGapResourceLimit("review_budget");
            ++out.kernel_calls;
            auto unrestricted=detail::AlignCentromereBlockDp(std::string_view(r.bases).substr(rp,rend-rp),std::string_view(q.bases).substr(qp,qend-qp),kernel);
            v.peak_kernel_bytes=std::max(v.peak_kernel_bytes,unrestricted.peak_kernel_bytes);v.dp_cells+=unrestricted.dp_cells;
            Position rr=rp,qq=qp;
            for(auto op:unrestricted.cigar){if(op.operation=='=')peer.push_back({rr,qq,op.length,0,0});if(op.operation!='I')rr+=op.length;if(op.operation!='D')qq+=op.length;}
          }
          Position rr=rp,qq=qp;
          for(auto op:gap.cigar){
            if(op.operation=='='){
              if(!opposing)confirmed.emplace_back(qq,qq+op.length);
              else for(const auto&other:peer){const auto first=std::max(qq,other.q),last=std::min(qq+op.length,other.q+other.n);if(first<last&&rr+first-qq==other.r+first-other.q)confirmed.emplace_back(first,last);}
            }
            Emit(v.new_cigar,op.operation,op.length);if(op.operation!='I')rr+=op.length;if(op.operation!='D')qq+=op.length;
          }
        }
        rp=rend;qp=qend;
      };
      for(const auto&s:consensus){solve(s.r,s.q);Emit(v.new_cigar,'=',s.n);confirmed.emplace_back(s.q,s.q+s.n);rp+=s.n;qp+=s.n;supported+=s.n;}
      solve(re,qe);
      JointFilterCoincidences(v,assembly,r,q,Union(std::move(confirmed)));
      // Confidence filtering may withhold pairs, but it cannot erase local
      // counterevidence and then present the surviving part as a new proof.
      // The filter spells an unconfirmed pair as equal I/D consumption.
      // That placeholder is not substitution evidence from the search.
      const bool substitution_conflict=composite && !core_review && !unplanned_review && JointNewSubstitutionCorner(v,assembly);
      if(constrained)JointMaskDiscordantComponents(v,assembly,current);
      if(composite)JointMaskOverlappingCopies(v,r,q,assembly,kernel,core_review || unplanned_review);
      if(core_review && !JointRestoreFlanks(v,current.final,r,q)){
        v.reason="existing_variable_context_not_preserved";
        v.seconds=Seconds(begin);out.reviews.push_back(std::move(v));continue;
      }
      for(const auto op:v.new_cigar)if(op.operation=='=')v.new_exact_bases+=op.length;
      v.new_score=detail::ScoreCentromereStitchCigar(std::string_view(r.bases).substr(rb,re-rb),std::string_view(q.bases).substr(qb,qe-qb),v.new_cigar,kernel.score_model);
      v.scored=true;v.engine="local_rare_chain_exact_gap_dp";
      const bool retained=constrained && JointRetainsExactPairs(assembly,v);
      const bool ordered=composite && JointOrderedContexts(r,q,v,kernel);
      const auto new_contexts=CrossPathContexts(r,rc,qc,v.new_cigar,rb,qb,v);
      v.old_unique_context_starts=prior_contexts.size();v.new_unique_context_starts=new_contexts.size();
      if(unplanned_review){
        if(!ordered || !std::includes(new_contexts.begin(),new_contexts.end(),prior_contexts.begin(),prior_contexts.end()))
          v.reason="unplanned_path_context_unresolved";
        else if(!JointSelectUnplannedRun(v,current,assembly,r,q,kernel))
          v.reason="unplanned_path_has_no_free_exact_run";
        else{
          const auto exact=v.query_end-v.query_begin;
          v.old_exact_bases=0;v.old_score=0;v.new_exact_bases=exact;
          v.new_score=detail::ScoreCentromereStitchCigar(
            std::string_view(r.bases).substr(v.reference_begin,exact),
            std::string_view(q.bases).substr(v.query_begin,exact),v.new_cigar,kernel.score_model);
          auto candidate=JointReplace(current,v,r,q);
          auto fixed=JointReplay(assembly,candidate,r,q,v);
          auto restitched=detail::AlignCentromereAnchoredPathsWithFixedGaps(r,q,candidate,copies,kernel,fixed);
          kernel.interruption_callback("joint-unplanned-publish");
          const auto actual=CrossFragments(restitched,r,q,v.reference_begin,v.reference_end,v.query_begin,v.query_end);
          if(!restitched.complete || restitched.aborted || actual.exact!=exact ||
             actual.mismatches || actual.score<=0)
            v.reason="unplanned_replay_unresolved";
          else{
            v.accepted=true;v.reason="strict_unplanned_run_gain";++out.accepted;
            out.anchors=std::move(candidate);out.stitched=std::move(restitched);
            used_r.emplace_back(v.reference_begin,v.reference_end);used_q.emplace_back(v.query_begin,v.query_end);
          }
        }
      }
      else if(substitution_conflict || (composite && !core_review && JointNewSubstitutionCorner(v,assembly)))v.reason="new_substitution_corner_unresolved";
      else if(constrained && !retained)v.reason="old_exact_pair_not_retained";
      else if(composite && !ordered)v.reason="ordered_context_competition_unresolved";
      else if(std::any_of(protected_copy_roots.begin(),protected_copy_roots.end(),[&](const auto&a){return !JointRetainsCopyRoot(a,v);}))v.reason="primary_multiplicity_requires_copy_review";
      else if(composite && !core_review && std::any_of(current.final.begin(),current.final.end(),[&](const auto&a){
        if(a.strand!=Strand::Forward||a.length<64||a.query_begin<qb||a.query_begin+a.length>qe||
           a.reference_begin<rb||a.reference_begin+a.length>re)return false;
        Position rr=rb,qq=qb;Length kept=0;
        for(auto op:v.new_cigar){
          if(op.operation=='='){const auto lo=std::max(qq,a.query_begin),hi=std::min(qq+op.length,a.query_begin+a.length);
            if(lo<hi&&rr+lo-qq==a.reference_begin+lo-a.query_begin)kept+=hi-lo;}
          if(op.operation!='I')rr+=op.length;
          if(op.operation!='D')qq+=op.length;
        }
        // A candidate that endorses the old locus cannot quietly move only
        // its remaining endpoint. Such edits belong to the constrained route.
        return kept && kept<a.length;
      }))v.reason="partial_old_placement_requires_endpoint_review";
      else if(composite && std::any_of(current.final.begin(),current.final.end(),[&](const auto&a){
        return a.strand==Strand::Forward && a.query_begin>=qb && a.query_begin+a.length<=qe &&
          a.reference_begin>=rb && a.reference_begin+a.length<=re &&
          JointVariableFlank(a,r,q) && !JointRetainsCopyRoot(a,v);
      }))v.reason="existing_variable_context_not_preserved";
      else if(!retained && JointQualityDominated(old,v.new_score,v.new_exact_bases,ShortDisagreementBurden(v.new_cigar)))v.reason="worse_fragment_consistency";
      else if(v.new_score<=v.old_score)v.reason=v.new_score==v.old_score?"equal_score_preserved":"lower_score_preserved";
      else if(v.new_exact_bases<=v.old_exact_bases)v.reason="no_exact_coverage_gain";
      else if(!std::includes(new_contexts.begin(),new_contexts.end(),prior_contexts.begin(),prior_contexts.end()))v.reason="existing_internal_context_not_preserved";
      else if(new_contexts.empty()&&consensus.size()<2)v.reason="insufficient_local_chain_context";
      else if(!exterior && !ordered && [&]{Length n=0;for(auto op:v.new_cigar)if(op.operation=='I'||op.operation=='D')n+=op.length;return n>=64;}())v.reason="no_independent_exterior_context";
      else if(!ordered && new_contexts.size()==prior_contexts.size() && JointLongGapWithoutContext(r,q,assembly,v))v.reason="long_repeat_gap_without_new_context";
      else if(new_contexts.empty() && JointMovableGap(r,q,v))v.reason="movable_repeat_gap_preserved";
      else{
        auto candidate=JointReplace(current,v,r,q);
        auto fixed=JointReplay(assembly,candidate,r,q,v);
        auto restitched=detail::AlignCentromereAnchoredPathsWithFixedGaps(r,q,candidate,copies,kernel,fixed);
        kernel.interruption_callback("joint-publish");
        if(!restitched.complete||restitched.aborted)v.reason="candidate_restitch_unresolved";
        else{
          const auto actual=CrossFragments(restitched,r,q,rb,re,qb,qe);
          if(actual.mismatches || actual.score<=old.score || actual.exact<=old.exact)v.reason="no_projected_score_exact_gain";
          else if(!retained && JointQualityDominated(old,actual.score,actual.exact,actual.short_disagreements))v.reason="projected_fragment_consistency";
          else{v.accepted=true;v.reason=core_review?"strict_internal_core_path_gain":retained?"strict_pair_retained_path_gain":ordered?"strict_ordered_context_path_gain":"strict_joint_context_path_gain";++out.accepted;out.anchors=std::move(candidate);out.stitched=std::move(restitched);used_r.emplace_back(rb,re);used_q.emplace_back(qb,qe);}
        }
      }
      (void)supported;
    }catch(const detail::StitchGapResourceLimit&error){v.reason=error.reason;v.peak_kernel_bytes=std::max(v.peak_kernel_bytes,error.peak_kernel_bytes);}
    catch(const std::bad_alloc&){v.reason="allocation_failure";}
    if(v.reason=="existing_internal_context_not_preserved"||v.reason=="existing_variable_context_not_preserved")
      marker_conflicts.push_back({rb,re,qb,qe});
    if(v.reason=="worse_fragment_consistency"||v.reason=="projected_fragment_consistency"){
      unresolved_contexts.push_back({rb,re,qb,qe});
      // Two disagreeing proposals in the same larger protected context must
      // not be resolved by taking only the favorable part of either path.
      // A single rejected proposal alone does not veto an alternative that
      // genuinely preserves all of that context's old locating markers.
      for(const auto&x:marker_conflicts)
        if(rb>=x.rb && re<=x.re && qb>=x.qb && qe<=x.qe)unresolved_contexts.push_back(x);
    }
    const bool accepted=v.accepted;
    v.seconds=Seconds(begin);out.reviews.push_back(std::move(v));
    if(composite && accepted){changed=true;break;}
    }
    // Accepted repairs can free reference positions that previously belonged
    // to short remote matches. Revisit the largest conflict first, retaining
    // all accepted footprints and the original shared resource accounting.
    if(composite && !changed){
      if(!constrained && !core_review && !unplanned_review){constrained=true;continue;}
      if(constrained && !marker_conflicts.empty()){constrained=false;core_review=true;continue;}
      if(!unplanned_review && !unplanned_contexts.empty()){
        constrained=false;core_review=false;unplanned_review=true;continue;
      }
      break;
    }
  }
  out.seconds=Seconds(started);return out;
#endif
}
} // namespace ramag
