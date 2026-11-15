#pragma once
#include "ramag/internal/centromere/direction.hpp"
#include <algorithm>
#include <limits>
#include <memory_resource>
#include <set>
#include <span>
#include <array>
#include <chrono>
#include <cmath>
#include <numeric>
#include <tuple>


// Related direction checks share a private module; their admission rules
// and template callback paths remain independent.

// orientation_connections
namespace ramag::detail {
// Endpoints are frozen direct-support intervals and original known blocks, not
// inferred gaps. Consecutive endpoints partition the remaining unknown space;
// thus a same-direction old block is a pivot, never overwritten or extrapolated.
template<class Emit,class Check>
void VisitConnectedOrientationGaps(std::span<const OrientationEndpoint> endpoints,
    std::span<const OrientationCoverage> labels,std::string_view query,Length limit,
    Emit emit,Check check,bool bridge_unknown=false,
    std::span<const CentromereQueryInterval> forward={},
    std::span<const CentromereQueryInterval> reverse={}) {
  for(std::size_t i=0;i<endpoints.size();++i) {
    check();const auto& right=endpoints[i];
    if(right.begin>=right.end || right.end>query.size() || (i && endpoints[i-1].end>right.begin))
      throw AlignmentError("invalid frozen orientation endpoints");
    if(!i)continue;
    const auto& left=endpoints[i-1];
    if(left.strand!=right.strand || left.end==right.begin)continue;
    CentromereOrientationGap gap;
    gap.begin=left.end;gap.end=right.begin;gap.strand=right.strand;
    gap.left_block_id=left.block_id;gap.right_block_id=right.block_id;
    gap.reason="gap_too_long";
    if(gap.end-gap.begin<=limit) {
      gap.reason="contains_N";
      if(query.substr(gap.begin,gap.end-gap.begin).find('N')==std::string_view::npos) {
        gap.reason="conflicting_or_uncovered_window";
        const unsigned bit=gap.strand==Strand::Forward?1U:2U;
        auto label=std::lower_bound(labels.begin(),labels.end(),gap.begin,
          [](const auto& x,auto p){return x.end<=p;});
        Position covered=gap.begin;
        for(;label!=labels.end() && covered<gap.end;++label) {
          // Require every position, including holes between windows, to have
          // exclusively same-direction evidence. Competition still vetoes.
          if(label->begin>covered || label->mask!=bit)break;
          covered=std::min(gap.end,label->end);
        }
        if(covered==gap.end){gap.accepted=true;gap.reason="continuous_window_gap";}
        else if(bridge_unknown) {
          // Keep the old accepted set. The fallback may cross missing coverage,
          // but a competing/opposite window anywhere in the whole gap vetoes it.
          auto candidate=std::lower_bound(labels.begin(),labels.end(),gap.begin,
            [](const auto& x,auto p){return x.end<=p;});
          bool veto=false;
          for(;candidate!=labels.end() && candidate->begin<gap.end;++candidate)
            if(candidate->mask & ~bit){veto=true;break;}
          if(!veto) {
            // Canonical query unions are sorted and disjoint; alternative copies
            // never add votes. Even one opposing matched base blocks this fallback.
            const auto opposite=gap.strand==Strand::Forward?reverse:forward;
            const auto match=std::lower_bound(opposite.begin(),opposite.end(),gap.begin,
              [](const auto& x,auto p){return x.end<=p;});
            gap.reason="opposite_exact_evidence";
            if(match==opposite.end() || match->begin>=gap.end) {
              gap.accepted=true;gap.reason="same_direction_unknown_gap";
            }
          }
        }
      }
    }
    emit(gap);
  }
}
} // namespace ramag::detail

// orientation_long_gap
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

// orientation_mixed_gap
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

// orientation_sparse_gap
namespace ramag::detail {
// Reuse the frozen old sparse graph and accepted exact vertices, not a second
// graph search. A recovered tip cannot make its next root reachable. Missing
// direction windows are allowed here, unlike the pure-window endpoint layers.
template<class Emit,class Check>
void VisitSparseOrientationGapEnds(
    std::span<const CentromereSparseOrientationGap> edges,
    std::span<const CentromereSparseOrientationInterval> accepted,
    std::span<const OrientationCoverage> labels,std::string_view query,Length limit,
    std::span<const CentromereQueryInterval> forward,
    std::span<const CentromereQueryInterval> reverse,
    std::span<const CentromereQueryInterval> occupied,Emit emit,Check check,
    bool prefer_stronger=false) {
  if(!limit)throw AlignmentError("sparse-gap endpoint recovery requires a positive limit");
  auto validate=[&](auto rows) {
    for(std::size_t i=0;i<rows.size();++i) {
      check();const auto& row=rows[i];
      if(row.begin>=row.end || row.end>query.size() || (i && rows[i-1].end>row.begin))
        throw AlignmentError("invalid sparse-gap interval index");
    }
  };
  validate(accepted);validate(labels);validate(forward);validate(reverse);validate(occupied);
  for(const auto& row:labels)if(row.mask>7U)throw AlignmentError("invalid sparse-gap window mask");
  for(const auto& edge:edges) {
    check();
    if(edge.begin>=edge.end || edge.end>query.size() || edge.parent_begin>edge.begin ||
       edge.parent_end<edge.end || edge.parent_end>query.size() ||
       (edge.strand!=Strand::Forward && edge.strand!=Strand::Reverse))
      throw AlignmentError("invalid frozen sparse-gap edge");
    // The original edge is tested before any occupied positions are subtracted.
    if(edge.accepted || edge.reason!="gap_too_long" || edge.end-edge.begin<=limit)continue;
    const Position begin=edge.begin,end=edge.end;
    const unsigned bit=edge.strand==Strand::Forward?1U:2U;
    const auto own=bit==1U?forward:reverse,opposite=bit==1U?reverse:forward;
    struct Tip {
      Position boundary{},root_begin{},root_end{};
      std::uint64_t source{};
      std::string_view reason{"no_frozen_sparse_root"};
      Length support{},competing_support{};
      bool scanned{},compared{};
    };
    auto qualify=[&](bool from_left) {
      const Position start=from_left?begin:end;
      Tip tip;tip.boundary=start;
      if((from_left && !begin) || (!from_left && end==query.size()))return tip;
      const Position pivot=from_left?begin-1:end;
      const auto root=std::lower_bound(own.begin(),own.end(),pivot,
        [](const auto& row,auto p){return row.end<=p;});
      if(root==own.end() || root->begin>pivot ||
         (from_left ? root->end!=begin : root->begin!=end))return tip;
      const auto source=std::lower_bound(accepted.begin(),accepted.end(),pivot,
        [](const auto& row,auto p){return row.end<=p;});
      if(source==accepted.end() || source->begin>pivot || source->inferred || !source->id ||
         source->strand!=edge.strand || source->parent_begin!=edge.parent_begin ||
         source->parent_end!=edge.parent_end || source->begin<root->begin || source->end>root->end)return tip;
      // A clipped accepted row is legal only at this ORIGINAL root boundary;
      // its clipped interior edge never becomes a new starting coordinate.
      tip.root_begin=root->begin;tip.root_end=root->end;tip.source=source->id;
      tip.reason="limit";return tip;
    };
    auto left=qualify(true),right=qualify(false);
    auto output=[&](bool from_left,const Tip& tip,bool accept,std::string_view reason) {
      CentromereSparseGapEnd row;
      row.begin=accept?(from_left?begin:tip.boundary):begin;
      row.end=accept?(from_left?tip.boundary:end):end;
      row.parent_begin=edge.parent_begin;row.parent_end=edge.parent_end;
      row.edge_begin=begin;row.edge_end=end;row.root_begin=tip.root_begin;row.root_end=tip.root_end;
      row.source_interval_id=tip.source;row.strand=edge.strand;row.from_left=from_left;row.accepted=accept;
      row.reason=from_left?"sparse_gap_endpoint_left_":"sparse_gap_endpoint_right_";
      row.reason+=reason;
      if(prefer_stronger) {
        row.meeting_compared=tip.compared;row.accepted_root_support=tip.support;
        row.competing_root_support=tip.competing_support;
        row.scan_begin=from_left?begin:tip.boundary;row.scan_end=from_left?tip.boundary:end;
        row.scan_reason=tip.scanned?tip.reason:"not_scanned";
      }
      emit(row);
    };
    // Do not convert a wholly supported long edge into a complete >limit bridge.
    auto cell=std::lower_bound(labels.begin(),labels.end(),begin,
      [](const auto& row,auto p){return row.end<=p;});
    Position covered=begin;
    for(;cell!=labels.end() && covered<end;++cell) {
      check();if(cell->begin>covered || cell->mask!=bit)break;
      covered=std::min(end,cell->end);
    }
    if(covered==end) {
      output(true,left,false,"whole_gap_pure_window");output(false,right,false,"whole_gap_pure_window");continue;
    }
    auto scan=[&](bool from_left,Tip& tip) {
      if(!tip.source)return;
      tip.scanned=true;
      const Position start=from_left?begin:end;
      const Position bound=from_left?begin+std::min(limit,end-begin):end-std::min(limit,end-begin);
      tip.boundary=bound;
      if(from_left) {
        auto label=std::lower_bound(labels.begin(),labels.end(),start,
          [](const auto& row,auto p){return row.end<=p;});
        for(;label!=labels.end() && label->begin<tip.boundary;++label) {
          check();
          if(label->mask & ~bit) {
            tip.boundary=std::max(start,label->begin);tip.reason="opposing_or_competitive_window";break;
          }
        }
        auto trim=[&](auto rows,std::string_view reason) {
          const auto row=std::lower_bound(rows.begin(),rows.end(),start,
            [](const auto& x,auto p){return x.end<=p;});
          if(row!=rows.end() && row->begin<tip.boundary) {
            tip.boundary=std::max(start,row->begin);tip.reason=reason;
          }
        };
        trim(opposite,"opposite_exact_evidence");trim(occupied,"occupied");
        const auto n=query.substr(start,tip.boundary-start).find_first_not_of("ACGT");
        if(n!=std::string_view::npos){tip.boundary=start+n;tip.reason="contains_N";}
      } else {
        auto label=std::lower_bound(labels.begin(),labels.end(),start,
          [](const auto& row,auto p){return row.begin<p;});
        while(label!=labels.begin()) {
          check();--label;if(label->end<=tip.boundary)break;
          if(label->mask & ~bit) {
            tip.boundary=std::min(start,label->end);tip.reason="opposing_or_competitive_window";break;
          }
        }
        auto trim=[&](auto rows,std::string_view reason) {
          auto row=std::lower_bound(rows.begin(),rows.end(),start,
            [](const auto& x,auto p){return x.begin<p;});
          if(row!=rows.begin() && (--row)->end>tip.boundary) {
            tip.boundary=std::min(start,row->end);tip.reason=reason;
          }
        };
        trim(opposite,"opposite_exact_evidence");trim(occupied,"occupied");
        const auto n=query.substr(tip.boundary,start-tip.boundary).find_last_not_of("ACGT");
        if(n!=std::string_view::npos){tip.boundary+=n+1;tip.reason="contains_N";}
      }
    };
    scan(true,left);scan(false,right);
    if(left.boundary>=right.boundary) {
      if(prefer_stronger && left.source && right.source && left.boundary>begin && right.boundary<end) {
        // accepted is disjoint and immutable. Visit only fragments intersecting
        // the complete root; clipped fragments and distinct accepted pieces all
        // contribute, while inferred rows and unrelated parents never do.
        auto support=[&](const Tip& tip) {
          Length total=0;
          auto row=std::lower_bound(accepted.begin(),accepted.end(),tip.root_begin,
            [](const auto& x,auto p){return x.end<=p;});
          for(;row!=accepted.end() && row->begin<tip.root_end;++row) {
            check();
            if(row->inferred || !row->id || row->strand!=edge.strand ||
               row->parent_begin!=edge.parent_begin || row->parent_end!=edge.parent_end)continue;
            total+=std::min(row->end,tip.root_end)-std::max(row->begin,tip.root_begin);
          }
          return total;
        };
        left.support=support(left);right.support=support(right);
        left.competing_support=right.support;right.competing_support=left.support;
        left.compared=true;right.compared=true;
        if(left.support==right.support) {
          output(true,left,false,"meeting_equal_support");output(false,right,false,"meeting_equal_support");
        } else {
          const bool keep_left=left.support>right.support;
          output(true,left,keep_left,keep_left?"meeting_kept_stronger_root":"meeting_weaker_root");
          output(false,right,!keep_left,keep_left?"meeting_weaker_root":"meeting_kept_stronger_root");
        }
        continue;
      }
      output(true,left,false,"ends_would_meet");output(false,right,false,"ends_would_meet");continue;
    }
    output(true,left,left.boundary>begin,left.reason);
    output(false,right,right.boundary<end,right.reason);
  }
}
} // namespace ramag::detail

// orientation_unified_ends
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


// Related direction checks share a private module; their admission rules
// and template callback paths remain independent.

// orientation_block_audit
namespace ramag::detail {
// Old chain envelopes remain provenance, not unconditional orientation truth.
// Audit only their internal, anchor-free long gaps. A pure opposite window cell
// retracts the old direction to unknown; it does not assign a new reference copy
// or turn the window's envelope into a new opposite-direction prediction.
// No paths or sequences are copied. Existing chains are reference-ascending;
// reverse chains can therefore be visited backwards in query order.
template<class EmitAudit,class EmitKnown,class Check>
void VisitKnownOrientationBlockAudit(std::span<const CentromereBlock> blocks,
    std::span<const OrientationCoverage> labels,Length query_length,Length gap_limit,
    EmitAudit audit,EmitKnown known,Check check,bool audit_all_supported_conflicts=false) {
  if(!gap_limit)throw AlignmentError("known-block audit requires a positive gap limit");
  for(std::size_t i=0;i<labels.size();++i) {
    check();const auto& row=labels[i];
    if(row.begin>=row.end || row.end>query_length || row.mask>7U ||
       (i && labels[i-1].end>row.begin))
      throw AlignmentError("invalid known-block audit coverage cells");
  }
  for(const auto& block:blocks) {
    check();if(!block.direction_known)continue;
    if(block.query_begin>block.query_end || block.query_end>query_length ||
       block.reference_begin>block.reference_end)
      throw AlignmentError("invalid known-block audit bounds");
    const auto& anchors=block.planning_anchors;
    for(std::size_t i=0;i<anchors.size();++i) {
      check();const auto& a=anchors[i];
      // Subtraction bounds the additions below and rejects overflow as well.
      if(!a.length || a.query_begin<block.query_begin || a.query_begin>block.query_end ||
         a.length>block.query_end-a.query_begin || a.reference_begin<block.reference_begin ||
         a.reference_begin>block.reference_end || a.length>block.reference_end-a.reference_begin)
        throw AlignmentError("planning anchor lies outside its audited block");
      if(i) {
        const auto& p=anchors[i-1];
        if(p.reference_begin+p.length>a.reference_begin ||
           (block.strand==Strand::Forward ? p.query_begin+p.length>a.query_begin :
             a.query_begin+a.length>p.query_begin))
          throw AlignmentError("planning anchors are not an ordered nonoverlapping chain");
      }
    }
    if(block.query_begin==block.query_end)continue;
    Position retained_begin=block.query_begin;
    const unsigned opposite=block.strand==Strand::Forward?2U:1U;
    auto anchor=[&](std::size_t i)->const CentromereAnchor& {
      return anchors[block.strand==Strand::Forward?i:anchors.size()-1-i];
    };
    for(std::size_t i=1;i<anchors.size();++i) {
      check();const auto& left=anchor(i-1);const auto& right=anchor(i);
      const Position gap_begin=left.query_begin+left.length,gap_end=right.query_begin;
      // The opt-in repair admits direct opposing-window evidence at any
      // positive gap size. Planning anchors themselves remain unchanged.
      if(gap_begin==gap_end || (!audit_all_supported_conflicts && gap_end-gap_begin<=gap_limit))continue;
      auto label=std::lower_bound(labels.begin(),labels.end(),gap_begin,
        [](const auto& row,auto p){return row.end<=p;});
      Position cut_begin=0,cut_end=0;bool pending=false;
      auto emit=[&] {
        if(!pending)return;
        if(retained_begin<cut_begin)
          known(CentromereOrientationKnownInterval{block.id,retained_begin,cut_begin,block.strand});
        audit(CentromereOrientationBlockAudit{block.id,cut_begin,cut_end,gap_begin,gap_end,block.strand});
        retained_begin=cut_end;pending=false;
      };
      for(;label!=labels.end() && label->begin<gap_end;++label) {
        check();if(label->mask!=opposite)continue;
        const auto begin=std::max(gap_begin,label->begin),end=std::min(gap_end,label->end);
        if(begin>=end)continue;
        if(pending && begin==cut_end)cut_end=end;
        else{emit();cut_begin=begin;cut_end=end;pending=true;}
      }
      emit();
    }
    if(retained_begin<block.query_end)
      known(CentromereOrientationKnownInterval{block.id,retained_begin,block.query_end,block.strand});
  }
}
} // namespace ramag::detail

// orientation_isolated_review
namespace ramag::detail {
// This does not remove a root or alter its rarity, chain, window or mapping.
// Only a bounded query-only gap can be emitted, from a frozen input snapshot.
template<class Nodes,class Neutral,class Emit,class EmitCounter,class Check>
void VisitIsolatedOppositeGaps(const Nodes& nodes,const Neutral& neutral,
    std::span<const CentromereUnifiedPivot> pivots,std::span<const OrientationCoverage> labels,
    std::string_view query,Length limit,Length neighborhood,
    std::span<const CentromereQueryInterval> known_forward,
    std::span<const CentromereQueryInterval> known_reverse,
    std::pmr::memory_resource* memory,Emit emit,EmitCounter emit_counter,Check check) {
  if(!limit || !neighborhood || !memory || neutral.size()!=nodes.size())
    throw AlignmentError("invalid isolated-opposite review configuration");
  auto intersects=[](auto rows,Position b,Position e) {
    const auto it=std::lower_bound(rows.begin(),rows.end(),b,[](const auto& x,auto p){return x.end<=p;});
    return it!=rows.end() && it->begin<e;
  };
  for(auto rows:{known_forward,known_reverse})for(std::size_t i=0;i<rows.size();++i) {
    check();if(rows[i].begin>=rows[i].end || rows[i].end>query.size() || (i && rows[i-1].end>=rows[i].begin))
      throw AlignmentError("isolated review requires merged known intervals");
  }
  for(const auto& x:known_forward)if(intersects(known_reverse,x.begin,x.end))
    throw AlignmentError("overlapping known directions");
  for(std::size_t i=0;i<labels.size();++i)if(labels[i].begin>=labels[i].end || labels[i].end>query.size() ||
      labels[i].mask>7U || (i && labels[i-1].end>labels[i].begin))throw AlignmentError("invalid isolated review labels");
  using Ids=std::pmr::vector<std::size_t>;
  std::array<Ids,2> roots{Ids(memory),Ids(memory)};
  std::array<std::pmr::vector<Position>,2> prefix{std::pmr::vector<Position>(memory),std::pmr::vector<Position>(memory)};
  for(std::size_t i=0;i<nodes.size();++i) {
    check();const auto& n=nodes[i];const auto& a=n.anchor;
    if(!a.length || a.query_begin>query.size() || a.length>query.size()-a.query_begin ||
       a.reference_begin>std::numeric_limits<Position>::max()-a.length || !std::isfinite(a.weight) || a.weight<0 ||
       (n.strand!=Strand::Forward && n.strand!=Strand::Reverse))throw AlignmentError("invalid normalized review root");
    if(!neutral[i])roots[n.strand==Strand::Forward?0U:1U].push_back(i);
  }
  for(unsigned s=0;s<2;++s) {
    auto& ids=roots[s];std::sort(ids.begin(),ids.end(),[&](auto a,auto b) {
      return std::pair{nodes[a].anchor.query_begin,a}<std::pair{nodes[b].anchor.query_begin,b};});
    Position e=0;for(auto i:ids){e=std::max(e,nodes[i].anchor.query_begin+nodes[i].anchor.length);prefix[s].push_back(e);}
  }
  auto visit=[&](unsigned sign,Position b,Position e,auto fn) {
    const auto& ids=roots[sign];const auto& ends=prefix[sign];
    auto k=static_cast<std::size_t>(std::upper_bound(ends.begin(),ends.end(),b)-ends.begin());
    for(;k<ids.size() && nodes[ids[k]].anchor.query_begin<e;++k) {
      check();const auto i=ids[k];const auto& a=nodes[i].anchor;
      if(a.query_begin+a.length>b && !fn(i))break;
    }
  };
  // Independent, two-sided collinearity, with the existing RaMA diagonal cost.
  auto coherent=[&](std::size_t i,std::size_t j) {
    if(i==j)return false;
    const auto* a=&nodes[i].anchor;const auto* b=&nodes[j].anchor;
    if(a->reference_begin>b->reference_begin)std::swap(a,b);
    const auto sign=nodes[i].strand;
    const Position au=sign==Strand::Forward?a->query_begin:query.size()-a->query_begin-a->length;
    const Position bu=sign==Strand::Forward?b->query_begin:query.size()-b->query_begin-b->length;
    if(a->reference_begin+a->length>b->reference_begin || au+a->length>bu)return false;
    if(b->reference_begin-a->reference_begin-a->length>neighborhood || bu-au-a->length>neighborhood)return false;
    const auto dr=b->reference_begin-a->reference_begin,dq=bu-au,delta=dr>dq?dr-dq:dq-dr;
    return !delta || 2.0*std::log2(static_cast<double>(delta)+1.0)<std::min(a->weight,b->weight);
  };
  using Pointers=std::pmr::vector<const CentromereUnifiedPivot*>;
  std::array<Pointers,4> ends{Pointers(memory),Pointers(memory),Pointers(memory),Pointers(memory)};
  auto position=[](const auto* p){return p->from_left?p->query_end:p->query_begin;};
  for(const auto& p:pivots) {
    check();if(!p.id || !p.source_id || p.source_kind>2 || p.query_begin>=p.query_end || p.query_end>query.size() ||
       p.reference_begin>=p.reference_end || p.reference_end-p.reference_begin!=p.query_end-p.query_begin ||
       (p.strand!=Strand::Forward && p.strand!=Strand::Reverse))throw AlignmentError("invalid isolated review pivot");
    ends[(p.strand==Strand::Forward?0U:2U)+(p.from_left?0U:1U)].push_back(&p);
  }
  for(auto& rows:ends)std::sort(rows.begin(),rows.end(),[&](auto a,auto b){return std::pair{position(a),a->id}<std::pair{position(b),b->id};});
  std::uint64_t next=0;
  for(unsigned sign=0;sign<2;++sign) {
    const auto known=sign==0?known_forward:known_reverse,opposite=sign==0?known_reverse:known_forward;
    const auto& lefts=ends[2*sign];const auto& rights=ends[2*sign+1];const unsigned forbidden=sign==0?6U:5U;
    for(std::size_t k=1;k<known.size();++k) {
      check();const auto& a=known[k-1];const auto& z=known[k];const auto b=a.end,e=z.begin;
      if(intersects(opposite,b,e))continue;
      CentromereLocalGapReview row;row.id=++next;row.begin=b;row.end=e;row.strand=sign==0?Strand::Forward:Strand::Reverse;
      auto decide=[&]() -> const char* {
        auto li=std::upper_bound(lefts.begin(),lefts.end(),b,[&](auto p,auto x){return p<position(x);});
        const auto ri=std::lower_bound(rights.begin(),rights.end(),e,[&](auto x,auto p){return position(x)<p;});
        if(li==lefts.begin() || ri==rights.end())return "missing_original_endpoint";
        --li;const auto l=position(*li),r=position(*ri);
        if(l<=a.begin || r>=z.end)return "missing_original_endpoint";
        li=std::lower_bound(lefts.begin(),lefts.end(),l,[&](auto x,auto p){return position(x)<p;});
        row.left_pivot=(*li)->id;row.right_pivot=(*ri)->id;row.corridor_begin=l;row.corridor_end=r;
        if(r-l>limit && r-l-limit>limit)return "original_endpoint_limit";
        if(intersects(opposite,l,r))return "opposite_known_direction";
        if(query.substr(l,r-l).find_first_not_of("ACGT")!=std::string_view::npos)return "contains_N";
        auto w=std::lower_bound(labels.begin(),labels.end(),l,[](const auto& x,auto p){return x.end<=p;});
        for(;w!=labels.end() && w->begin<r;++w)if(w->mask & forbidden)return "opposing_or_competitive_window";
        bool in_gap=false,clear=true;
        visit(1U-sign,l,r,[&](auto i) {
          const auto& root=nodes[i].anchor;
          in_gap=in_gap || (root.query_begin<e && root.query_begin+root.length>b);
          CentromereOppositeReview counter;counter.review_id=row.id;counter.normalized_id=i+1;
          counter.reference_begin=root.reference_begin;counter.query_begin=root.query_begin;
          counter.length=root.length;counter.weight=root.weight;
          if(root.query_begin<l || root.query_begin+root.length>r)counter.verdict=1;
          else if(intersects(opposite,root.query_begin,root.query_begin+root.length))counter.verdict=2;
          else if(root.length>=(*li)->query_end-(*li)->query_begin || root.length>=(*ri)->query_end-(*ri)->query_begin)counter.verdict=3;
          else {
            auto low=root.query_begin-std::min(neighborhood,root.query_begin);
            const auto end=root.query_begin+root.length;
            auto high=end+std::min(neighborhood,static_cast<Length>(query.size())-end);
            if(low)--low;
            if(high<query.size())++high; // Include exact neighborhood-boundary partners.
            visit(1U-sign,low,high,[&](auto j) {
              if(coherent(i,j)){counter.neighbor_id=j+1;counter.verdict=4;return false;}return true;
            });
          }
          emit_counter(counter);clear=clear && !counter.verdict;return true;
        });
        if(!in_gap)return "no_opposite_root_in_gap";
        if(!clear)return "opposite_evidence_retains_veto";
        row.accepted=true;return "isolated_opposite_locally_outsupported";
      };
      row.reason=decide();emit(row);
    }
  }
}
} // namespace ramag::detail

// orientation_known_counterchain
namespace ramag::detail {
// A final, query-only subtraction. No search callback, inferred pivot or
// re-planning path is exposed: every decision reads the frozen upstream result.
void ReviewKnownBlockCounterchains(std::string_view query,
    CentromereBlockResult& out,const CentromereBlockOptions& options);
} // namespace ramag::detail

// orientation_local_direction_repair
namespace ramag::detail {
// One final query-direction view. Candidate roots, challenged gaps, protections
// and unknown space are frozen before any decision; repairs never feed back.
void RepairInheritedDirections(std::string_view query,
    CentromereBlockResult& out,const CentromereBlockOptions& options);
} // namespace ramag::detail
