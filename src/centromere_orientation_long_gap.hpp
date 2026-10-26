#pragma once
#include "centromere_orientation_connections.hpp"

namespace ramag::detail {
// A separate, one-pass overlay on the original long gaps. The spans passed here
// are frozen before this pass: accepted ends can never become new pivots. The
// short-gap specialization uses the same scan but restricts the ORIGINAL gap
// to <=limit and clips to existing occupancy, never to inferred new endpoints. Known
// directions are post-audit; exact support is the nonneutral root-query union.
template<bool ShortGaps=false,class Emit,class Check>
void VisitLongOrientationGapEnds(std::span<const OrientationEndpoint> endpoints,
    std::span<const OrientationCoverage> labels,std::string_view query,Length limit,
    std::span<const CentromereQueryInterval> known_forward,
    std::span<const CentromereQueryInterval> known_reverse,
    std::span<const CentromereQueryInterval> forward,
    std::span<const CentromereQueryInterval> reverse,
    std::span<const CentromereQueryInterval> occupied,Emit emit,Check check) {
  if(!limit)throw AlignmentError("long-gap endpoint recovery requires a positive limit");
  auto validate=[&](auto rows) {
    for(std::size_t i=0;i<rows.size();++i) {
      check();const auto& row=rows[i];
      if(row.begin>=row.end || row.end>query.size() || (i && rows[i-1].end>row.begin))
        throw AlignmentError("invalid long-gap interval index");
    }
  };
  validate(endpoints);validate(labels);validate(known_forward);validate(known_reverse);
  validate(forward);validate(reverse);validate(occupied);
  for(const auto& row:labels)if(row.mask>7U)throw AlignmentError("invalid long-gap window mask");
  auto contains=[](auto rows,Position p) {
    const auto it=std::lower_bound(rows.begin(),rows.end(),p,
      [](const auto& row,auto x){return row.end<=x;});
    return it!=rows.end() && it->begin<=p;
  };
  for(std::size_t i=1;i<endpoints.size();++i) {
    check();const auto& left=endpoints[i-1];const auto& right=endpoints[i];
    const auto length=right.begin-left.end;
    if(left.strand!=right.strand || !length || (ShortGaps ? length>limit : length<=limit))continue;
    if(left.strand!=Strand::Forward && left.strand!=Strand::Reverse)
      throw AlignmentError("invalid long-gap endpoint strand");
    const Position begin=left.end,end=right.begin;
    const unsigned bit=left.strand==Strand::Forward?1U:2U;
    const auto own=bit==1U?forward:reverse,opposite=bit==1U?reverse:forward;
    const auto known=bit==1U?known_forward:known_reverse;
    auto output=[&](bool from_left,Position b,Position e,bool accepted,std::string_view reason) {
      CentromereOrientationGap row;
      row.begin=accepted?b:begin;row.end=accepted?e:end;row.strand=left.strand;
      row.left_block_id=left.block_id;row.right_block_id=right.block_id;
      row.accepted=accepted;
      if constexpr(ShortGaps)row.reason=from_left?"short_gap_endpoint_left_":"short_gap_endpoint_right_";
      else row.reason=from_left?"long_gap_endpoint_left_":"long_gap_endpoint_right_";
      row.reason+=reason;emit(row);
    };
    auto reject=[&](std::string_view reason) {
      output(true,begin,end,false,reason);output(false,begin,end,false,reason);
    };
    // Audit-created cut edges are not endpoints. Both original boundary bases
    // must still have their old direction and actual root exact-match support.
    if(!contains(known,left.end-1) || !contains(own,left.end-1)) {
      reject("left_pivot_not_supported");continue;
    }
    if(!contains(known,right.begin) || !contains(own,right.begin)) {
      reject("right_pivot_not_supported");continue;
    }
    // A wholly covered long gap is a different policy from terminal recovery.
    // Do not silently turn two short extensions into a >limit complete bridge.
    auto cell=std::lower_bound(labels.begin(),labels.end(),begin,
      [](const auto& row,auto p){return row.end<=p;});
    Position covered=begin;
    for(;cell!=labels.end() && covered<end;++cell) {
      check();if(cell->begin>covered || cell->mask!=bit)break;
      covered=std::min(end,cell->end);
    }
    if constexpr(!ShortGaps) {
      if(covered==end){reject("whole_gap_pure_window");continue;}
    }
    // Short gaps may have pure windows but an internal opposite exact match.
    // Scan only from the original ends; that match remains an uncrossable veto.
    struct Tip {Position boundary{};std::string_view reason{"limit"};};
    auto scan=[&](bool from_left) {
      const Position start=from_left?begin:end;
      const Position bound=from_left?begin+std::min(limit,end-begin):end-std::min(limit,end-begin);
      Tip tip{start,"limit"};
      if(from_left) {
        auto it=std::lower_bound(labels.begin(),labels.end(),start,
          [](const auto& row,auto p){return row.end<=p;});
        while(tip.boundary<bound) {
          check();
          if(it==labels.end() || it->begin>tip.boundary || it->mask!=bit) {
            tip.reason=it!=labels.end() && it->begin<=tip.boundary && it->mask?
              "opposing_or_competitive_window":"uncovered_window";break;
          }
          tip.boundary=std::min(bound,it->end);++it;
        }
        auto trim=[&](auto rows,std::string_view reason) {
          const auto match=std::lower_bound(rows.begin(),rows.end(),start,
            [](const auto& row,auto p){return row.end<=p;});
          if(match!=rows.end() && match->begin<tip.boundary) {
            tip.boundary=std::max(start,match->begin);tip.reason=reason;
          }
        };
        trim(opposite,"opposite_exact_evidence");trim(occupied,"occupied");
        const auto n=query.substr(start,tip.boundary-start).find('N');
        if(n!=std::string_view::npos){tip.boundary=start+n;tip.reason="contains_N";}
      } else {
        auto it=std::lower_bound(labels.begin(),labels.end(),start,
          [](const auto& row,auto p){return row.begin<p;});
        while(tip.boundary>bound) {
          check();
          if(it==labels.begin()) {tip.reason="uncovered_window";break;}
          --it;
          if(it->end<tip.boundary || it->mask!=bit) {
            tip.reason=it->end>=tip.boundary && it->mask?
              "opposing_or_competitive_window":"uncovered_window";break;
          }
          tip.boundary=std::max(bound,it->begin);
        }
        auto trim=[&](auto rows,std::string_view reason) {
          auto match=std::lower_bound(rows.begin(),rows.end(),start,
            [](const auto& row,auto p){return row.begin<p;});
          if(match!=rows.begin() && (--match)->end>tip.boundary) {
            tip.boundary=std::min(start,match->end);tip.reason=reason;
          }
        };
        trim(opposite,"opposite_exact_evidence");trim(occupied,"occupied");
        const auto n=query.substr(tip.boundary,start-tip.boundary).rfind('N');
        if(n!=std::string_view::npos){tip.boundary+=n+1;tip.reason="contains_N";}
      }
      return tip;
    };
    const auto l=scan(true),r=scan(false);
    if(l.boundary>=r.boundary){reject("ends_would_meet");continue;}
    output(true,begin,l.boundary,l.boundary>begin,l.reason);
    output(false,r.boundary,end,r.boundary<end,r.reason);
  }
}
} // namespace ramag::detail
