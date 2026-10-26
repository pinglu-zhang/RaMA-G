#pragma once
#include "centromere_orientation_unified_ends.hpp"
#include <array>

namespace ramag::detail {
// Read-only snapshot: an accepted review never supplies another endpoint or
// changes a window label. Nearest facing ORIGINAL pivots dominate farther ones
// for both distance and vetoes, so no endpoint-pair enumeration is necessary.
template<class Emit,class Check>
void VisitLocalCompetitiveGaps(std::span<const CentromereUnifiedPivot> pivots,
    std::span<const OrientationCoverage> labels,std::string_view query,Length limit,
    std::span<const CentromereQueryInterval> forward,std::span<const CentromereQueryInterval> reverse,
    std::span<const CentromereQueryInterval> known_forward,std::span<const CentromereQueryInterval> known_reverse,
    std::pmr::memory_resource* memory,Emit emit,Check check) {
  if(!limit || !memory)throw AlignmentError("local gap review requires positive limit and memory resource");
  auto validate=[&](auto rows) {
    for(std::size_t i=0;i<rows.size();++i) {
      check();const auto& row=rows[i];
      if(row.begin>=row.end || row.end>query.size() || (i && rows[i-1].end>=row.begin))
        throw AlignmentError("local review requires merged interval indices");
    }
  };
  validate(forward);validate(reverse);validate(known_forward);validate(known_reverse);
  for(std::size_t i=0;i<labels.size();++i) {
    const auto& row=labels[i];
    if(row.begin>=row.end || row.end>query.size() || row.mask>7U || (i && labels[i-1].end>row.begin))
      throw AlignmentError("invalid local review window coverage");
  }
  auto intersects=[](auto rows,Position begin,Position end) {
    const auto row=std::lower_bound(rows.begin(),rows.end(),begin,[](const auto& x,auto p){return x.end<=p;});
    return row!=rows.end() && row->begin<end;
  };
  for(const auto& row:known_forward)if(intersects(known_reverse,row.begin,row.end))
    throw AlignmentError("opposite known direction overlap");
  using Index=std::pmr::vector<const CentromereUnifiedPivot*>;
  std::array<Index,4> index{Index(memory),Index(memory),Index(memory),Index(memory)};
  for(const auto& p:pivots) {
    check();
    if(!p.id || !p.source_id || p.source_kind>2 || p.query_begin>=p.query_end || p.query_end>query.size() ||
       p.reference_begin>=p.reference_end || p.reference_end-p.reference_begin!=p.query_end-p.query_begin ||
       (p.strand!=Strand::Forward && p.strand!=Strand::Reverse))
      throw AlignmentError("invalid original local-review endpoint");
    index[(p.strand==Strand::Forward?0U:2U)+(p.from_left?0U:1U)].push_back(&p);
  }
  auto position=[](const auto* p){return p->from_left?p->query_end:p->query_begin;};
  for(auto& ids:index)std::sort(ids.begin(),ids.end(),[&](auto a,auto b){
    return std::pair{position(a),a->id}<std::pair{position(b),b->id};});
  std::uint64_t id=0;
  for(unsigned sign=0;sign<2;++sign) {
    const auto known=sign==0?known_forward:known_reverse;
    const auto opposite_known=sign==0?known_reverse:known_forward,opposite_exact=sign==0?reverse:forward;
    const unsigned own_bit=1U<<sign,opposite_bit=sign==0?2U:1U;
    const auto& lefts=index[2*sign];const auto& rights=index[2*sign+1];
    for(std::size_t i=1;i<known.size();++i) {
      check();const auto& left=known[i-1];const auto& right=known[i];
      const Position begin=left.end,end=right.begin;
      // Adjacent known components must be separated only by unknown positions.
      if(intersects(opposite_known,begin,end))continue;
      CentromereLocalGapReview row;row.id=++id;row.begin=begin;row.end=end;
      row.strand=sign==0?Strand::Forward:Strand::Reverse;
      auto decide=[&]() -> const char* {
        auto li=std::upper_bound(lefts.begin(),lefts.end(),begin,[&](auto p,auto a){return p<position(a);});
        const auto ri=std::lower_bound(rights.begin(),rights.end(),end,[&](auto a,auto p){return position(a)<p;});
        if(li==lefts.begin() || ri==rights.end())return "missing_original_endpoint";
        --li;const Position lb=position(*li),rb=position(*ri);
        // End base must belong to the immediately adjacent same-sign component.
        if(lb<=left.begin || rb>=right.end)return "missing_original_endpoint";
        // Deterministic smallest ID when duplicate positions are supplied.
        li=std::lower_bound(lefts.begin(),lefts.end(),lb,[&](auto a,auto p){return position(a)<p;});
        row.left_pivot=(*li)->id;row.right_pivot=(*ri)->id;
        row.corridor_begin=lb;row.corridor_end=rb;
        const Length distance=rb-lb;
        if(distance>limit && distance-limit>limit)return "original_endpoint_limit";
        if(intersects(opposite_exact,lb,rb))return "opposite_exact_inside_corridor";
        if(intersects(opposite_known,lb,rb))return "opposite_known_inside_corridor";
        if(query.substr(lb,rb-lb).find_first_not_of("ACGT")!=std::string_view::npos)return "contains_N";
        bool competitive=false;Position covered=begin;
        auto w=std::lower_bound(labels.begin(),labels.end(),lb,[](const auto& x,auto p){return x.end<=p;});
        for(;w!=labels.end() && w->begin<rb;++w) {
          check();
          if(w->mask & opposite_bit)return "opposite_advantage_window";
          if(w->end<=begin || w->begin>=end)continue;
          competitive=competitive || (w->mask & 4U);
          if(w->mask & own_bit) {
            if(w->begin>covered)return "same_advantage_not_continuous";
            covered=std::max(covered,std::min(end,w->end));
          }
        }
        if(!competitive)return "not_competitive_gap";
        if(covered!=end)return "same_advantage_not_continuous";
        row.accepted=true;return "local_corridor_clear";
      };
      row.reason=decide();emit(row);
    }
  }
}
} // namespace ramag::detail
