#pragma once
#include "centromere_orientation_connections.hpp"
#include <limits>
#include <memory_resource>
#include <set>

namespace ramag::detail {
inline bool QualifyUnifiedPlanningEnd(const CentromereAnchor& anchor,std::uint64_t block_id,
    Strand strand,bool from_left,std::span<const CentromereOrientationKnownInterval> known) {
  if(!block_id || !anchor.length || anchor.query_begin>std::numeric_limits<Position>::max()-anchor.length)return false;
  const Position p=from_left?anchor.query_begin+anchor.length-1:anchor.query_begin;
  const auto row=std::lower_bound(known.begin(),known.end(),p,[](const auto& x,auto v){return x.end<=v;});
  return row!=known.end() && row->begin<=p && row->block_id==block_id && row->strand==strand;
}
template<class Accepted>
bool QualifyUnifiedRootEnd(const CentromereAnchor& anchor,Strand strand,bool from_left,bool neutral,
    std::span<const CentromereQueryInterval> own,std::span<const Accepted> accepted) {
  if(neutral || !anchor.length || anchor.query_begin>std::numeric_limits<Position>::max()-anchor.length)return false;
  const Position begin=anchor.query_begin,end=begin+anchor.length,p=from_left?end-1:begin;
  const auto root=std::lower_bound(own.begin(),own.end(),p,[](const auto& x,auto v){return x.end<=v;});
  if(root==own.end() || root->begin>p)return false;
  const auto row=std::lower_bound(accepted.begin(),accepted.end(),p,[](const auto& x,auto v){return x.end<=v;});
  return row!=accepted.end() && row->begin<=p && !row->inferred && row->id && row->strand==strand;
}

// All inputs are immutable. Scan distance is measured from the ORIGINAL real
// endpoint even while traversing an already accepted same-direction prefix.
// Final subtraction is separate, and outputs never become new pivots.
template<class EmitScan,class EmitInterval,class Check>
void VisitUnifiedOrientationEnds(std::span<const CentromereUnifiedPivot> pivots,
    std::span<const OrientationCoverage> labels,std::string_view query,Length limit,
    std::span<const CentromereQueryInterval> forward,std::span<const CentromereQueryInterval> reverse,
    std::span<const CentromereQueryInterval> known_forward,std::span<const CentromereQueryInterval> known_reverse,
    std::span<const CentromereQueryInterval> occupied,std::pmr::memory_resource* memory,
    EmitScan emit_scan,EmitInterval emit_interval,Check check) {
  if(!limit || !memory)throw AlignmentError("unified endpoints require positive limit and memory resource");
  auto validate=[&](auto rows) {
    for(std::size_t i=0;i<rows.size();++i) {
      check();const auto& row=rows[i];
      if(row.begin>=row.end || row.end>query.size() || (i && rows[i-1].end>row.begin))
        throw AlignmentError("invalid unified endpoint interval index");
    }
  };
  validate(labels);validate(forward);validate(reverse);validate(known_forward);validate(known_reverse);validate(occupied);
  for(const auto& row:labels)if(row.mask>7U)throw AlignmentError("invalid unified endpoint window mask");
  struct Event {Position p;std::uint64_t id;unsigned sign;bool enter;};
  std::pmr::vector<Event> events(memory);
  std::uint64_t next_scan=0;
  for(const auto& pivot:pivots) {
    check();
    if(!pivot.id || !pivot.source_id || pivot.source_kind>2 || pivot.query_begin>=pivot.query_end ||
       pivot.query_end>query.size() || pivot.reference_begin>=pivot.reference_end ||
       pivot.reference_end-pivot.reference_begin!=pivot.query_end-pivot.query_begin ||
       (pivot.strand!=Strand::Forward && pivot.strand!=Strand::Reverse))
      throw AlignmentError("invalid unified original endpoint");
    const bool right=pivot.from_left;
    const Position start=right?pivot.query_end:pivot.query_begin;
    Position boundary=right?start+std::min(limit,static_cast<Length>(query.size())-start):start-std::min(limit,start);
    const unsigned sign=pivot.strand==Strand::Forward?0U:1U,bit=1U<<sign;
    const auto opposite=sign==0?reverse:forward,opposite_known=sign==0?known_reverse:known_forward;
    std::string_view reason="limit";
    if(right) {
      auto label=std::lower_bound(labels.begin(),labels.end(),start,[](const auto& x,auto p){return x.end<=p;});
      for(;label!=labels.end() && label->begin<boundary;++label) {
        check();if(label->mask & ~bit){boundary=std::max(start,label->begin);reason="opposing_or_competitive_window";break;}
      }
      auto trim=[&](auto rows,std::string_view why) {
        const auto row=std::lower_bound(rows.begin(),rows.end(),start,[](const auto& x,auto p){return x.end<=p;});
        if(row!=rows.end()) {
          const Position clipped=std::max(start,row->begin);
          if(clipped<boundary){boundary=clipped;reason=why;}
        }
      };
      trim(opposite,"opposite_exact_evidence");trim(opposite_known,"opposite_known_direction");
      const auto n=query.substr(start,boundary-start).find_first_not_of("ACGT");
      if(n!=std::string_view::npos){boundary=start+n;reason="contains_N";}
    } else {
      auto label=std::lower_bound(labels.begin(),labels.end(),start,[](const auto& x,auto p){return x.begin<p;});
      while(label!=labels.begin()) {
        check();--label;if(label->end<=boundary)break;
        if(label->mask & ~bit){boundary=std::min(start,label->end);reason="opposing_or_competitive_window";break;}
      }
      auto trim=[&](auto rows,std::string_view why) {
        auto row=std::lower_bound(rows.begin(),rows.end(),start,[](const auto& x,auto p){return x.begin<p;});
        if(row!=rows.begin()) {
          const Position clipped=std::min(start,(--row)->end);
          if(clipped>boundary){boundary=clipped;reason=why;}
        }
      };
      trim(opposite,"opposite_exact_evidence");trim(opposite_known,"opposite_known_direction");
      const auto n=query.substr(boundary,start-boundary).find_last_not_of("ACGT");
      if(n!=std::string_view::npos){boundary+=n+1;reason="contains_N";}
    }
    CentromereUnifiedScan scan{++next_scan,pivot.id,std::min(start,boundary),std::max(start,boundary),pivot.strand,std::string(reason)};
    emit_scan(scan);
    if(scan.begin<scan.end){events.push_back({scan.begin,scan.id,sign,true});events.push_back({scan.end,scan.id,sign,false});}
  }
  for(const auto& row:occupied){events.push_back({row.begin,0,2,true});events.push_back({row.end,0,2,false});}
  std::sort(events.begin(),events.end(),[](const auto& a,const auto& b){return a.p<b.p;});
  std::pmr::set<std::uint64_t> plus(memory),minus(memory);
  std::int64_t used=0;Position previous=0;std::uint64_t next_interval=0;
  for(std::size_t i=0;i<events.size();) {
    check();const Position p=events[i].p;
    if(previous<p && !used && (plus.empty()!=minus.empty())) {
      const bool positive=!plus.empty();const auto& active=positive?plus:minus;
      emit_interval(CentromereUnifiedInterval{++next_interval,previous,p,
        positive?Strand::Forward:Strand::Reverse,*active.begin()});
    }
    do {
      const auto& e=events[i];
      if(e.sign==2)used+=e.enter?1:-1;
      else {auto& active=e.sign==0?plus:minus;if(e.enter)active.insert(e.id);else active.erase(e.id);}
      ++i;
    } while(i<events.size() && events[i].p==p);
    previous=p;
  }
}
} // namespace ramag::detail
