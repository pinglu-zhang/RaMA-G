#pragma once
#include <array>
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace ramag {
enum class RepairWork : std::size_t {
    Tasks, CandidateCigar, LocalBases, SeedHits, ChainEdges, DpWork,
    ReplayColumns, PlacementBases, Count
};
enum class RepairLimitReason : std::uint8_t { None, Work, WallTime, Watchdog, Memory };
inline constexpr std::size_t RepairWorkCount=static_cast<std::size_t>(RepairWork::Count);
inline constexpr std::string_view RepairWorkName(RepairWork kind) noexcept {
    constexpr std::array<std::string_view,RepairWorkCount> names{
        "tasks","candidate_cigar","local_bases","seed_hits","chain_edges",
        "dp_work","replay_columns","placement_bases"};
    const auto i=static_cast<std::size_t>(kind);
    return i<names.size()?names[i]:"invalid";
}
inline constexpr std::string_view RepairLimitName(RepairLimitReason reason) noexcept {
    switch(reason) {
        case RepairLimitReason::None:return "none";
        case RepairLimitReason::Work:return "work";
        case RepairLimitReason::WallTime:return "wall_time";
        case RepairLimitReason::Watchdog:return "watchdog";
        case RepairLimitReason::Memory:return "memory";
    }
    return "invalid";
}
struct RepairBudgetLimits {
    std::array<std::uint64_t,RepairWorkCount> work{
        2048,100000000,50000000,20000000,50000000,500000000,
        100000000,std::numeric_limits<std::uint64_t>::max()};
    std::uint64_t auxiliary_bytes{512ULL<<20};
    std::uint64_t diagnostic_bytes{8ULL<<20};
    std::uint64_t diagnostic_rows{10000};
    double task_seconds{10.0};
};
struct RepairBudgetStatistics {
    std::array<std::uint64_t,RepairWorkCount> used{},limits{};
    std::uint64_t memory_current_bytes{},memory_peak_bytes{},memory_limit_bytes{};
    std::uint64_t diagnostic_bytes{},diagnostics_kept{},diagnostics_dropped{};
    std::uint64_t reference_scan_passes{};
    double seconds{},wall_limit_seconds{};
    bool exhausted{};
    RepairLimitReason reason{RepairLimitReason::None};
    std::string phase;
};
class RepairBudgetExceeded : public std::runtime_error {
public:
    RepairBudgetExceeded(RepairLimitReason reason,std::string_view phase)
        :std::runtime_error(std::string("repair budget exhausted: ")+
          std::string(RepairLimitName(reason))+" during "+std::string(phase)),reason_(reason){}
    [[nodiscard]] RepairLimitReason reason()const noexcept{return reason_;}
private:
    RepairLimitReason reason_;
};
class RepairBudget;
class RepairMemoryLease {
public:
    RepairMemoryLease()=default;
    RepairMemoryLease(const RepairMemoryLease&)=delete;
    RepairMemoryLease& operator=(const RepairMemoryLease&)=delete;
    RepairMemoryLease(RepairMemoryLease&& other)noexcept
        :owner_(std::exchange(other.owner_,nullptr)),bytes_(std::exchange(other.bytes_,0)){}
    RepairMemoryLease& operator=(RepairMemoryLease&& other)noexcept {
        if(this!=&other){Reset();owner_=std::exchange(other.owner_,nullptr);bytes_=std::exchange(other.bytes_,0);}return *this;
    }
    ~RepairMemoryLease(){Reset();}
    void Grow(std::uint64_t total_bytes);
    void Reset()noexcept;
    [[nodiscard]] std::uint64_t bytes()const noexcept{return bytes_;}
private:
    friend class RepairBudget;
    explicit RepairMemoryLease(RepairBudget* owner):owner_(owner){}
    RepairBudget* owner_{};
    std::uint64_t bytes_{};
};
class RepairTaskGuard {
public:
    RepairTaskGuard()=default;
    void Check()const;
    void Check(RepairWork work,std::uint64_t amount)const;
private:
    friend class RepairBudget;
    using Clock=std::chrono::steady_clock;
    RepairTaskGuard(RepairBudget* owner,Clock::time_point deadline,std::string_view phase)
        :owner_(owner),deadline_(deadline),phase_(phase){}
    RepairBudget* owner_{};
    Clock::time_point deadline_{};
    std::string phase_;
};
class RepairBudget {
public:
    using Clock=std::chrono::steady_clock;
    explicit RepairBudget(double wall_seconds,RepairBudgetLimits limits={})
        :limits_(std::move(limits)),start_(Clock::now()),wall_seconds_(wall_seconds) {
        if(!(wall_seconds>0.0)||wall_seconds>300.0||
           !(limits_.task_seconds>0.0)||limits_.task_seconds>10.0)
            throw std::invalid_argument("invalid cumulative repair time budget");
        // Public hosts may tighten a quota; they cannot bypass this model's
        // frozen cumulative caps by supplying a custom ledger.
        const RepairBudgetLimits maximum{};
        for(std::size_t i=0;i<RepairWorkCount;++i)
            if(limits_.work[i]>maximum.work[i])
                throw std::invalid_argument("repair work limit exceeds frozen cap");
        if(limits_.auxiliary_bytes>maximum.auxiliary_bytes||
           limits_.diagnostic_bytes>maximum.diagnostic_bytes||
           limits_.diagnostic_rows>maximum.diagnostic_rows)
            throw std::invalid_argument("repair memory or diagnostic limit exceeds frozen cap");
        stats_.limits=limits_.work;stats_.memory_limit_bytes=limits_.auxiliary_bytes;
        stats_.wall_limit_seconds=wall_seconds_;
    }
    RepairBudget(const RepairBudget&)=delete;
    RepairBudget& operator=(const RepairBudget&)=delete;
    void Check(std::string_view phase={}) {
        std::lock_guard lock(mutex_);CheckLocked(phase);
    }
    void Consume(RepairWork kind,std::uint64_t amount) {
        std::lock_guard lock(mutex_);CheckLocked(RepairWorkName(kind));
        const auto i=static_cast<std::size_t>(kind);
        if(i>=RepairWorkCount)throw std::invalid_argument("invalid repair work kind");
        if(amount>limits_.work[i]-stats_.used[i])FailLocked(RepairLimitReason::Work,RepairWorkName(kind));
        stats_.used[i]+=amount;
    }
    [[nodiscard]] RepairTaskGuard Task(std::string_view phase) {
        Consume(RepairWork::Tasks,1);return Watch(phase);
    }
    [[nodiscard]] RepairTaskGuard Watch(std::string_view phase) {
        Check(phase);
        return RepairTaskGuard(this,Clock::now()+
            std::chrono::duration_cast<Clock::duration>(std::chrono::duration<double>(limits_.task_seconds)),phase);
    }
    [[nodiscard]] RepairMemoryLease Memory(std::uint64_t bytes) {
        RepairMemoryLease result(this);result.Grow(bytes);return result;
    }
    void ClaimReferenceScan(std::uint64_t bases) {
        std::lock_guard lock(mutex_);CheckLocked("reference_scan");
        if(stats_.reference_scan_passes)FailLocked(RepairLimitReason::Work,"reference_scan_passes");
        const auto i=static_cast<std::size_t>(RepairWork::PlacementBases);
        if(bases>limits_.work[i]-stats_.used[i])FailLocked(RepairLimitReason::Work,"placement_bases");
        stats_.used[i]+=bases;++stats_.reference_scan_passes;
    }
    // Diagnostics are deliberately dropped rather than exhausting scientific work.
    // Their conservative storage reservation is held until this invocation ends.
    bool KeepDiagnostic(std::uint64_t bytes) {
        std::lock_guard lock(mutex_);
        if(stats_.diagnostics_kept>=limits_.diagnostic_rows||
           bytes>limits_.diagnostic_bytes-stats_.diagnostic_bytes||
           bytes>limits_.auxiliary_bytes-stats_.memory_current_bytes){
            ++stats_.diagnostics_dropped;return false;
        }
        stats_.diagnostic_bytes+=bytes;++stats_.diagnostics_kept;
        stats_.memory_current_bytes+=bytes;
        stats_.memory_peak_bytes=std::max(stats_.memory_peak_bytes,stats_.memory_current_bytes);
        return true;
    }
    [[nodiscard]] RepairBudgetStatistics Snapshot()const {
        std::lock_guard lock(mutex_);auto out=stats_;
        out.seconds=std::chrono::duration<double>(Clock::now()-start_).count();return out;
    }
private:
    friend class RepairMemoryLease;
    friend class RepairTaskGuard;
    [[noreturn]] void FailLocked(RepairLimitReason reason,std::string_view phase) {
        if(!stats_.exhausted){stats_.exhausted=true;stats_.reason=reason;stats_.phase=phase;}
        throw RepairBudgetExceeded(stats_.reason,stats_.phase);
    }
    void CheckLocked(std::string_view phase) {
        if(stats_.exhausted)throw RepairBudgetExceeded(stats_.reason,stats_.phase);
        if(std::chrono::duration<double>(Clock::now()-start_).count()>=wall_seconds_)
            FailLocked(RepairLimitReason::WallTime,phase);
    }
    void TaskCheck(Clock::time_point deadline,std::string_view phase) {
        std::lock_guard lock(mutex_);CheckLocked(phase);
        if(Clock::now()>=deadline)FailLocked(RepairLimitReason::Watchdog,phase);
    }
    void Acquire(std::uint64_t bytes) {
        std::lock_guard lock(mutex_);CheckLocked("auxiliary_memory");
        if(bytes>limits_.auxiliary_bytes-stats_.memory_current_bytes)
            FailLocked(RepairLimitReason::Memory,"auxiliary_memory");
        stats_.memory_current_bytes+=bytes;
        stats_.memory_peak_bytes=std::max(stats_.memory_peak_bytes,stats_.memory_current_bytes);
    }
    void Release(std::uint64_t bytes)noexcept {
        std::lock_guard lock(mutex_);stats_.memory_current_bytes-=bytes;
    }
    RepairBudgetLimits limits_;
    Clock::time_point start_;
    double wall_seconds_;
    mutable std::mutex mutex_;
    RepairBudgetStatistics stats_;
};
inline void RepairMemoryLease::Grow(std::uint64_t total_bytes) {
    if(total_bytes<=bytes_)return;
    if(!owner_)throw std::logic_error("inactive repair memory lease");
    owner_->Acquire(total_bytes-bytes_);bytes_=total_bytes;
}
inline void RepairMemoryLease::Reset()noexcept {
    if(owner_){owner_->Release(bytes_);owner_=nullptr;bytes_=0;}
}
inline void RepairTaskGuard::Check()const {
    if(owner_)owner_->TaskCheck(deadline_,phase_);
}
inline void RepairTaskGuard::Check(RepairWork work,std::uint64_t amount)const {
    Check();if(owner_)owner_->Consume(work,amount);
}
} // namespace ramag
