#pragma once
#include "centromere_orientation_connections.hpp"

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
