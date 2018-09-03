/*
 * A simple multithreaded interval-based bnb solver
 */

/*
 * File:   tutorialbnb.cpp
 * Author: yamchenko.y.v
 *
 * Created on January 3, 2018, 5:10 PM
 */

#include <iostream>
#include <limits>
#include <random>
#include <algorithm>
#include <vector>
#include <iterator>
#include <functional>
#include <thread>
#include <chrono>
#include <forward_list>
#include <mutex>
#include <memory>
#include <atomic>
#include <numeric>
#include <unordered_set>
#include <condition_variable>

#include "logger.hpp"

#include <common/parbench.hpp>

using BM = Benchmark<double>;
using Box = std::vector<Interval<double>>;

constexpr char gKnownRecord[] = "knrec";

static double gEps = 0.1;

static bool gKnrec = false;

static double gSubsSplitCoeff = 0.5;

static double gStepsSplitCoeff = 0.5;

static int gMtSubsLimit = 10;

static int gMaxStepsTotal = 1000000;

static int gProcs = 4;

static int gMtStepsLimit = 1000;

std::atomic<double> gRecord;

Logger logger("paral_solver.log");

/**
 * BnB state
 */
struct ThreadHolder {

    enum Status {
        IS_READY, IS_PROCESSING, IS_FINISHED
    };

    ThreadHolder(ThreadPool* tPool, Status stat = IS_READY) : mThreadPool(tPool), mStatus(stat) {
    }

    ThreadHolder(const ThreadHolder& st) : mRecordVal(st.mRecordVal),
    mRecord(st.mRecord),
    mPool(st.mPool),
    mRemainingSteps(st.mRemainingSteps),
    mStatus(st.mStatus.load()),
    mMutex() {}

    ThreadHolder& operator=(const ThreadHolder& st) {
        mRecordVal = st.mRecordVal;
        mRecord = st.mRecord;
        mPool.assign(st.mPool.begin(), st.mPool.end());
        mRemainingSteps = st.mRemainingSteps;
        mStatus.store(st.mStatus.load());
        return *this;
    }

    void merge(const std::shared_ptr<ThreadHolder> s) {
        mRemainingSteps += s->mRemainingSteps;
        s->mRemainingSteps = 0;
        if (s->mRecordVal < mRecordVal) {
            mRecordVal = s->mRecordVal;
            mRecord = s->mRecord;
        }
        mPool.insert(mPool.end(), s->mPool.begin(), s->mPool.end());
    }

    bool hasResources() {
        return ! mPool.empty() && (mRemainingSteps > 0);
    }

    void printResourses() {
        std::cout << "Remaining task count: " << mPool.size() << "\n"
                "Remainng steps count: " << mRemainingSteps << "\n";
    }

    void assignTaskTo(std::shared_ptr<ThreadHolder> st) {
        *st = *this;
        mRemainingSteps = 0;
        mPool.clear();
    }

    std::shared_ptr<ThreadHolder> trySplit(std::shared_ptr<ThreadHolder> s = nullptr) {
        std::lock_guard<std::mutex> lock(mMutex);
        const int pool_size = mPool.size();
        if (pool_size <= gMtSubsLimit || mRemainingSteps <= gMtStepsLimit) {
            return nullptr;
        }
        if (s == nullptr) {
            s = std::make_shared<ThreadHolder>();
        }
        const int mv_step_count = mRemainingSteps * gStepsSplitCoeff;
        mRemainingSteps -= mv_step_count;
        s->mRemainingSteps = mv_step_count;
        s->mRecord = mRecord;
        s->mRecordVal = mRecordVal;
        const int mv_sub_count = pool_size * gSubsSplitCoeff;
        auto mv_end_iter = mPool.begin() + mv_sub_count;
        s->mPool.assign(mPool.begin(), mv_end_iter);
        mPool.erase(mPool.begin(), mv_end_iter);
        return s;
    }

    Box getSub() {
        std::lock_guard<std::mutex> lock(mMutex);
        const Box box = mPool.back();
        mPool.pop_back();
        mRemainingSteps--;
        return box;
    }

    void clear() {
        mRecordVal = std::numeric_limits<double>::max();
        mRecord.clear();
        mPool.clear();
        mRemainingSteps = 0;
    }

    void setReady() {
        clear();
        mStatus = IS_READY;
    }

    void setProcessing() {
        mStatus = IS_PROCESSING;
    }

    void setFinished() {
        mStatus = IS_FINISHED;
    }

    bool isFinished() const {
        return (this->mStatus == IS_FINISHED) ? true : false;
    }

    bool isProcessing() const {
        return (this->mStatus == IS_PROCESSING) ? true : false;
    }

    bool isReady() const {
        return (this->mStatus == IS_READY) ? true : false;
    }

    void run(ThreadPool* const tp) {

    }

    double mRecordVal = 0.0;
    std::vector<double> mRecord;
    std::vector<Box> mPool;
    int mRemainingSteps = 0;
    std::atomic<Status> mStatus;
    std::mutex mMutex;
    ThreadPool* mThreadPool;
};

class Notifier {
public:

    Notifier() : mNotificationCount(0) {
    }

    void notify() {
        mCV.notify_one();
        std::atomic_fetch_add(&mNotificationCount, 1);
    }

    void resolve() {
        int nc = mNotificationCount.load();
        if (nc > 0) {
            std::atomic_fetch_add(&mNotificationCount, 1);
        }
    }

    void wait_notification(std::shared_ptr<ThreadHolder> s) {
        std::unique_lock<std::mutex> lock(s->mMutex);
        mCV.wait_for(lock, mWaitPeriod, [&]() {
            return mNotificationCount.load() != 0;
        });
    }

private:
    std::condition_variable mCV;
    std::atomic<int> mNotificationCount;
    std::chrono::duration<int64_t, std::milli> mWaitPeriod = std::chrono::milliseconds(100);
};

Notifier gNotifier;

class ThreadPool final {
public:
//    ThreadPool(const size_t pool_size, const BM& bm, bool knrec = false) :
//        mReadyStatesIndexes(pool_size), mList(pool_size), mBenchmark(bm) {
//        std::iota(mReadyStatesIndexes.begin(), mReadyStatesIndexes.end(), 0);

//        init(knrec);
//    }

    ThreadPool(const size_t pool_size, const BM& bm, bool knrec = false) :
            mReadyStatesIndexes(pool_size), mList(pool_size), mBenchmark(bm) {
            std::iota(mReadyStatesIndexes.begin(), mReadyStatesIndexes.end(), 0);

            init(knrec);
        }

    size_t popReadyThreadHolderIdx() {
        if (this->mReadyStatesIndexes.empty()) {
            return -1;
        }
        return mReadyStatesIndexes.pop_back();
    }

    void addReadyThreadHolderIdx(const size_t idx) {
        mReadyStatesIndexes.push_back(idx);
    }

    size_t activeCount() const {
        return mActiveThreadCount;
    }

    ///TODO: подумать
    void loadThread() {
        ++mActiveThreadCount;
    }

    void releaseThread(size_t idx) {
        this->mReadyStatesIndexes.push_back(idx);
        --mActiveThreadCount;
    }

    void run() {
        size_t idx = popReadyThreadHolderIdx();
        this->mList[idx].run(this);
        this->loadThread();
    }

private:
    void init(bool knrec) {
        const int dim = this->mBenchmark.getDim();
        Box ibox;

        /// TODO: Вынести в Box
        for (size_t i = 0; i < dim; i++) {
            ibox.emplace_back(bm.getBounds()[i].first, bm.getBounds()[i].second);
        }

        ThreadHolder& s = mList[0];
        s->mPool.push_back(ibox);
        s->mRecordVal = std::numeric_limits<double>::max();
        s->mRemainingSteps = gMaxStepsTotal;

        if (knrec) {
            s->mRecordVal = bm.getGlobMinY();
        } else {
            s->mRecordVal = std::numeric_limits<double>::max();
        }
    }

    std::vector<ThreadHolder> mList;
    std::vector<size_t> mReadyStatesIndexes;
    size_t mActiveThreadCount = 0;
    BM& mBenchmark;
};

ThreadPool* gThreadPool;

std::ostream& operator<<(std::ostream & out, const std::shared_ptr<ThreadHolder> s) {
    out << "\"recval\" : " << s->mRecordVal << "\n";
    out << "\"record\" : [";
    for (int i = 0; i < s->mRecord.size(); i++) {
        out << s->mRecord[i];
        if (i != s->mRecord.size() - 1)
            out << ", ";
    }
    out << "]\n";
    out << "\"max steps\" :" << s->mRemainingSteps << "\n";
    return out;
}

double len(const Interval<double>& I) {
    return I.rb() - I.lb();
}

void split(const Box& ibox, std::shared_ptr<ThreadHolder> s) {
    auto result = std::max_element(ibox.begin(), ibox.end(),
            [](const Interval<double>& f, const Interval<double>& s) {
                return len(f) < len(s);
            });
    const int i = result - ibox.begin();
    const double maxlen = len(ibox[i]);
    Box b1(ibox);
    Interval<double> ilow(ibox[i].lb(), ibox[i].lb() + 0.5 * maxlen);
    b1[i] = ilow;
    Box b2(ibox);
    Interval<double> iupper(ibox[i].lb() + 0.5 * maxlen, ibox[i].rb());
    b2[i] = iupper;

    // LOCK AT EVERY ITERATION !!!
    std::lock_guard<std::mutex> lock(s->mMutex);
    s->mPool.push_back(std::move(b1));
    s->mPool.push_back(std::move(b2));
}

void getCenter(const Box& ibox, std::vector<double>& c) {
    const int n = ibox.size();
    for (int i = 0; i < n; i++) {
        c[i] = 0.5 * (ibox[i].lb() + ibox[i].rb());
    }
}

void solveSerial(std::shared_ptr<ThreadHolder> s, const BM& bm) {
    const int dim = bm.getDim();
    std::vector<double> c(dim);

    while (s->hasResources()) {
        Box b = s->getSub();
        getCenter(b, c);
        double v = bm.calcFunc(c);
        double rv = gRecord.load();
        while (v < rv) {
            if (gRecord.compare_exchange_strong(rv, v)) {
                s->mRecordVal = v;
                s->mRecord = c;
            }
        }
        auto lb = bm.calcInterval(b).lb();
        if (lb <= gRecord.load() - gEps) {
            split(b, s);
        }
    }
    s->setFinished();
    gNotifier.notify();
}

void runThread(std::shared_ptr<ThreadHolder> s, const BM& bm) {
    gThreadPool.mActiveThreadCount++;
    s->setProcessing();
    std::thread init_thread(solveSerial, s, std::ref(bm));
    init_thread.detach();
}

bool try_assign_run(std::shared_ptr<ThreadHolder> sender, std::shared_ptr<ThreadHolder> receiver, const BM& bm) {
    if (sender->hasResources() && sender->mRemainingSteps > gMtStepsLimit) {
        sender->assignTaskTo(receiver);
        runThread(receiver, bm);
        return true;
    } else {
        return false;
    }
}

void solve(std::shared_ptr<ThreadHolder> init_s, const BM& bm) {
    std::shared_ptr<ThreadHolder> first_s = std::make_shared<ThreadHolder>();
    init_s->assignTaskTo(first_s);
    gThreadPool.PushFront(first_s);

    runThread(first_s, bm);


    std::cout << "Remaining steps: " << init_s->mRemainingSteps << std::endl;
    if (init_s->hasResources()) {
        logger << "Has resources" << std::endl;
        init_s->printResourses();
        solveSerial(init_s, bm);
    }

    gThreadPool.clear();

    std::cout << "Problem is solved" << std::endl;
}

void initTask(const BM& bm) {
    gRecord.store(std::numeric_limits<double>::max());
    gThreadPool = new ThreadPool(gProcs, bm, gKnrec);
}

double findMin(const BM& bm) {
    initTask(bm);

    std::chrono::time_point<std::chrono::steady_clock> start, end;
    start = std::chrono::steady_clock::now();
#if 0
    solveSerial(&s, bm);
#else
    solve(s, bm);
#endif
    end = std::chrono::steady_clock::now();
    int mseconds = (std::chrono::duration_cast<std::chrono::microseconds> (end - start)).count();
    std::cout << "Time: " << mseconds << " microsecond\n";
    const int steps = gMaxStepsTotal - s->mRemainingSteps;
    std::cout << "Time per subproblem: " << ((mseconds > 0) ? ((double) mseconds / (double) steps) : 0) << " miscroseconds." << std::endl;
    if (s->mRemainingSteps == 0) {
        std::cout << "Failed to converge in " << gMaxStepsTotal << " steps\n";
    } else {
        std::cout << "Converged in " << steps << " steps\n";
    }

    std::cout << "BnB found = " << gRecord.load() << std::endl;
    std::cout << " at x [ ";
    std::copy(s->mRecord.begin(), s->mRecord.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << "]\n";
    return gRecord.load();

    delete gThreadPool;
}

void printInfo() {
    std::cout << "Record is lock-free: " << gRecord.is_lock_free() << std::endl;
    if (gKnrec == gKnownRecord) {
        std::cout << "Init global record: true" << std::endl;
    } else {
        std::cout << "Init global record: false" << std::endl;
    }
    std::cout << "Eps: " << gEps << std::endl;
    std::cout << "Max steps count: " << gMaxStepsTotal << std::endl;
    std::cout << "Max thread count: " << gProcs << std::endl;
    std::cout << "Available thread count: " << std::thread::hardware_concurrency() << std::endl;
    std::cout << "Steps split coeffitient count: " << gStepsSplitCoeff << std::endl;
    std::cout << "Subs split coeffitient count: " << gSubsSplitCoeff << std::endl;
    std::cout << "Steps limit: " << gMtStepsLimit << std::endl;
    std::cout << "Subs limit: " << gMtSubsLimit << std::endl << std::endl;
}

bool testBench(const BM& bm) {
    std::cout << "*************Testing benchmark**********" << std::endl;
    printInfo();
    std::cout << bm;

    bool res_flag = true;
    double v = findMin(bm);
    double diff = v - bm.getGlobMinY();
    if (diff > gEps) {
        res_flag = false;
        std::cout << "BnB failed for " << bm.getDesc() << " benchmark " << std::endl;
    }
    std::cout << "the difference is " << v - bm.getGlobMinY() << std::endl;
    std::cout << "****************************************" << std::endl << std::endl;

    return res_flag;
}

void printHelp(std::string bin_name) {
    std::cout << "Usage: " << bin_name << " '<name_of_bench>' <eps max_steps> <virtual_procs_number> <split_steps_limit> "
                                          "<split_subs_limit> <split_steps_coeff> <split_subs_coeff>\n\n";
    std::cout << "to run all tests run\n";
    std::cout << bin_name << std::endl << std::endl;
    std::cout << "to list benchmarks run:\n";
    std::cout << bin_name << " list\n\n";
    std::cout << "to see this message run\n";
    std::cout << bin_name << " --help\n";
}

main(int argc, char** argv) {
    logger << "Start program\n";
    std::string bench;
    ParBenchmarks<double> tests;
    if ((argc == 2) && (std::string(argv[1]) == std::string("list"))) {
        for (auto b : tests) {
            std::cout << b->getDesc() << "\n";
        }
        return 0;
    } else if (argc == 10) {
        bench = argv[1];
        std::string knrec_param = argv[2];
        if (knrec_param == gKnownRecord) {
            gKnrec = true;
        } else {
            gKnrec = false;
        }

        gEps = atof(argv[3]);
        gMaxStepsTotal = atoi(argv[4]);
        gProcs = atoi(argv[5]);
        gStepsSplitCoeff = atof(argv[6]);
        gSubsSplitCoeff = atof(argv[7]);
        gMtStepsLimit = atoi(argv[8]);
        gMtSubsLimit = atoi(argv[9]);
    } else {
        std::cerr << "Usage: " << argv[0] << " name_of_bench knrec|unknrec eps max_steps thread_count steps_split_coeff subs_split_coeff step_limit subs_limit\n";
        std::cerr << "or to list benchmarks run:\n";
        std::cerr << argv[0] << " list\n";
        return -1;
    }
    for (auto bm : tests) {
        if (bench == bm->getDesc())
            testBench(*bm);
    }
    logger << "Program end\n" << std::endl;
}
