#pragma once
#include "centromere_orientation_connections.hpp"

namespace ramag::detail {
// An independent overlay on original, opposite-direction endpoints. Each tip
// keeps its own direction; neither the other tip nor inferred output is a pivot.
template<class Emit,class Check>
void VisitMixedOrientationGapEnds(std::span<const OrientationEndpoint> endpoints,
    std::span<const OrientationCoverage> labels,std::string_view query,Length limit,
    std::span<const CentromereQueryInterval> known_forward,
    std::span<const CentromereQueryInterval> known_reverse,
    std::span<const CentromereQueryInterval> forward,
    std::span<const CentromereQueryInterval> reverse,
    std::span<const CentromereQueryInterval> occupied,Emit emit,Check check) {
  if(!limit)throw AlignmentError("mixed-gap endpoint recovery requires a positive limit");
  auto validate=[&](auto rows) {
    for(std::size_t i=0;i<rows.size();++i) {
      check();const auto& row=rows[i];
      if(row.begin>=row.end || row.end>query.size() || (i && rows[i-1].end>row.begin))
        throw AlignmentError("invalid mixed-gap interval index");
    }
  };
  validate(endpoints);validate(labels);validate(known_forward);validate(known_reverse);
  validate(forward);validate(reverse);validate(occupied);
  for(const auto& row:endpoints)
    if(row.strand!=Strand::Forward && row.strand!=Strand::Reverse)
      throw AlignmentError("invalid mixed-gap endpoint strand");
  for(const auto& row:labels)if(row.mask>7U)throw AlignmentError("invalid mixed-gap window mask");
  auto contains=[](auto rows,Position p) {
    const auto it=std::lower_bound(rows.begin(),rows.end(),p,
      [](const auto& row,auto x){return row.end<=x;});
    return it!=rows.end() && it->begin<=p;
  };
  for(std::size_t i=1;i<endpoints.size();++i) {
    check();const auto& left=endpoints[i-1];const auto& right=endpoints[i];
    const Position begin=left.end,end=right.begin;
    if(left.strand==right.strand || begin==end)continue;
    struct Tip {Position boundary{};std::string_view reason{"limit"};};
    auto scan=[&](bool from_left) {
      const auto strand=from_left?left.strand:right.strand;
      const unsigned bit=strand==Strand::Forward?1U:2U;
      const auto own=bit==1U?forward:reverse,opposite=bit==1U?reverse:forward;
      const auto known=bit==1U?known_forward:known_reverse;
      const Position start=from_left?begin:end,pivot=from_left?begin-1:end;
      Tip tip{start,"limit"};
      if(!contains(known,pivot) || !contains(own,pivot)) {
        tip.reason="pivot_not_supported";return tip;
      }
      const Position bound=from_left?begin+std::min(limit,end-begin):end-std::min(limit,end-begin);
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
        const auto n=query.substr(start,tip.boundary-start).find_first_not_of("ACGT");
        if(n!=std::string_view::npos){tip.boundary=start+n;tip.reason="contains_N";}
      } else {
        auto it=std::lower_bound(labels.begin(),labels.end(),start,
          [](const auto& row,auto p){return row.begin<p;});
        while(tip.boundary>bound) {
          check();
          if(it==labels.begin()){tip.reason="uncovered_window";break;}
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
        const auto n=query.substr(tip.boundary,start-tip.boundary).find_last_not_of("ACGT");
        if(n!=std::string_view::npos){tip.boundary+=n+1;tip.reason="contains_N";}
      }
      return tip;
    };
    auto output=[&](bool from_left,Position b,Position e,bool accepted,std::string_view reason) {
      CentromereOrientationGap row;
      row.begin=accepted?b:begin;row.end=accepted?e:end;
      row.strand=from_left?left.strand:right.strand;
      row.left_block_id=left.block_id;row.right_block_id=right.block_id;row.accepted=accepted;
      row.reason=from_left?"mixed_gap_endpoint_left_":"mixed_gap_endpoint_right_";
      row.reason+=reason;emit(row);
    };
    const auto l=scan(true),r=scan(false);
    // Do not manufacture a direction boundary by letting execution order decide
    // which tip wins. Even touching tips leave the whole proposal unresolved.
    if(l.boundary>=r.boundary) {
      output(true,begin,end,false,"ends_would_meet");
      output(false,begin,end,false,"ends_would_meet");continue;
    }
    output(true,begin,l.boundary,l.boundary>begin,l.reason);
    output(false,r.boundary,end,r.boundary<end,r.reason);
  }
}
} // namespace ramag::detail
