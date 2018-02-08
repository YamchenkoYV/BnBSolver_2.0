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

#include <testfuncs/benchmarks.hpp>

using BM = Benchmark<double>;
using Box = std::vector<Interval<double>>;
struct State;

const static int gProcs = 2;

const static int gMtStepsLimit = 1000;
const static int gMtSubsLimit = 10;

const static int gMaxStepsTotal = 10000000;

struct ThreadList {

    State* getReadyState();

   std::forward_list<State*>& operator()() { return mList; }
   size_t Size() const { return mThreadListSize; }
   size_t ActiveCount() const { return mActiveThreadCount; }

   void PushFront(State* st) {
       ++mThreadListSize;
       mList.push_front(st);
       std::cout << "ThreadList size "<< mThreadListSize << std::endl;
   }

   void clear();

    std::forward_list< State* > mList;

    size_t mThreadListSize = 0;

    size_t mActiveThreadCount = 0;
};

ThreadList gThreadList;

struct State {
    enum Status { IS_READY, IS_PROCESSING, IS_FINISHED };

    State(Status stat = IS_READY) : mStatus(stat) {}
    State(const State& st) : mRecordVal(st.mRecordVal),
        mRecord(st.mRecord),
        mPool(st.mPool),
        mMaxSteps(st.mMaxSteps),
        mSteps(st.mSteps),
        mStatus(st.mStatus),
        mutex() {}

    State& operator=(const State& st) {
        mRecordVal = st.mRecordVal;
        mRecord = st.mRecord;
        mPool.assign(st.mPool.begin(), st.mPool.end());
        mMaxSteps = st.mMaxSteps;
        mSteps = st.mSteps;
        mStatus = st.mStatus;

        return *this;
    }

    void merge(const State& s) {
        mMaxSteps += s.mMaxSteps - s.mSteps;
        mSteps += s.mSteps;
        if (s.mRecordVal < mRecordVal) {
            mRecordVal = s.mRecordVal;
            mRecord = s.mRecord;
        }
        mPool.insert(mPool.end(), s.mPool.begin(), s.mPool.end());
    }

    bool hasResources() { return !mPool.empty() && mMaxSteps != 0; }

    void assignTaskTo(State* st) {
        *st = *this;
        clear();
    }

    State* try_split(double q, State* s = nullptr) {
        std::lock_guard<std::mutex> lock(mutex);

        const int remMaxSteps = mMaxSteps - mSteps;

        if (mPool.size() <= gMtSubsLimit || remMaxSteps <= gMtStepsLimit) {
            return nullptr;
        }

        if (s == nullptr) {
            s = new State();
        }

        mMaxSteps = mSteps + remMaxSteps / 2;
        s->mMaxSteps = remMaxSteps - mMaxSteps;

        s->mRecord = mRecord;
        s->mRecordVal = mRecordVal;

        const int mvSubCount = gMtSubsLimit * q;
        auto mvSubBegIter = mPool.end() - mvSubCount;

        s->mPool.assign(mPool.begin(), mvSubBegIter);
        mPool.erase(mPool.begin(), mvSubBegIter);

        return s;
    }

    Box getSub() {
        //std::lock_guard<std::mutex> lock(mutex);
        const Box box = mPool.back();
        mPool.pop_back();

        return box;
    }

    void clear() {
        mRecordVal = std::numeric_limits<double>::max();
        mRecord.clear();
        mPool.clear();
        mMaxSteps = 0;
        mSteps = 0;
    }

    void setReady() {
        clear();
        mStatus = IS_READY;
    }

    void setProcessing() { mStatus = IS_PROCESSING; }
    void setFinished() { mStatus = IS_FINISHED; }

    bool isFinished() const { return (this->mStatus == IS_FINISHED) ? true : false; }
    bool isProcessing() const { return (this->mStatus == IS_PROCESSING) ? true : false; }
    bool isReady() const { return (this->mStatus == IS_READY) ? true : false; }


    double mRecordVal;

    std::vector<double> mRecord;

    std::vector<Box> mPool;

    int mMaxSteps;

    int mSteps = 0;

    Status mStatus;

    std::mutex mutex;
};

State* ThreadList::getReadyState() {
    for (State* st : mList) {
        if (st->isReady()) {
            return st;
        }
    }

    return nullptr;
}

void ThreadList::clear() {
    for(auto it = mList.begin(); it != mList.end(); ++it) {
        State* cur_state = *it;
        delete(cur_state);
    }
}

std::ostream& operator<<(std::ostream & out, const State s) {
    out << "\"recval\" : " << s.mRecordVal << "\n";
    out << "\"record\" : [";
    for (int i = 0; i < s.mRecord.size(); i++) {
        out << s.mRecord[i];
        if (i != s.mRecord.size() - 1)
            out << ", ";
    }
    out << "]\n";
    out << "\"steps\" :" << s.mSteps << "\n";
    out << "\"max steps\" :" << s.mMaxSteps << "\n";
    return out;
}

double len(const Interval<double>& I) {
    return I.rb() - I.lb();
}

void split(const Box& ibox, std::vector<Box>& v) {
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
    v.push_back(std::move(b1));
    v.push_back(std::move(b2));
}

void getCenter(const Box& ibox, std::vector<double>& c) {
    const int n = ibox.size();
    for (int i = 0; i < n; i++) {
        c[i] = 0.5 * (ibox[i].lb() + ibox[i].rb());
    }
}

void solveSerial(State* s, const BM& bm, double eps) {

    const int dim = bm.getDim();
    std::vector<double> c(dim);
    while (!s->mPool.empty()) {
        std::lock_guard<std::mutex> lock(s->mutex);
        s->mSteps ++;
        Box b = s->getSub(); //Lock state's mutex
        getCenter(b, c);
        double v = bm.calcFunc(c);
        if (v < s->mRecordVal) {
            s->mRecordVal = v;
            s->mRecord = c;
        }
        auto lb = bm.calcInterval(b).lb();
        if (lb <= s->mRecordVal - eps) {
            split(b, s->mPool);
        }
        if (s->mSteps >= s->mMaxSteps)
            break;
    }

    s->setFinished();
}

bool isProblemSolved(){ return gThreadList.ActiveCount() == 0; }

void runThread(State* s, const BM& bm, double eps) {
    gThreadList.mActiveThreadCount = 1;
    s->setProcessing();
    std::thread init_thread(solveSerial, s, std::ref(bm), eps);
    init_thread.detach();
}

void solve(State& init_s, const BM& bm, double eps, double b_coeff) {
    if (b_coeff >= 1) {
        std::cerr << "Warn: balance coefficient should be less than 1."
        << " b_coeff set to 0.5" << std::endl;
        b_coeff = 0.5;
    }

    State* first_s = new State();
    init_s.assignTaskTo(first_s);
    gThreadList.PushFront(first_s);

    runThread(first_s, bm, eps);
    
    while(! isProblemSolved()){
        for( auto iter = gThreadList().begin() ; iter != gThreadList().end(); ++iter)
        {
            State* cur_state = *iter;
            if (cur_state->isReady()) {
                /*if (! init_s.mPool.empty() && init_s.hasResources()) {
                    init_s.assignTaskTo(cur_state);
                    runThread(cur_state, bm, eps);
                }*/
                continue;
            }
            
            if( cur_state->isFinished() )
            {
                init_s.merge(*cur_state);
                cur_state->setReady();
                --gThreadList.mActiveThreadCount;
            } else if( gThreadList.mActiveThreadCount < gProcs)
            {
                //Load balancing procedures
                    State* new_state;
                    if (gThreadList.mActiveThreadCount == gThreadList.Size()) {
                        new_state = cur_state->try_split(b_coeff);
                        if( new_state == nullptr )
                            continue;

                        gThreadList.PushFront(new_state);
                    } else {
                        new_state = gThreadList.getReadyState();
                        if( new_state == nullptr )
                            continue;

                        new_state = cur_state->try_split(b_coeff, new_state);
                        if( new_state == nullptr )
                            continue;
                    }

                    runThread(new_state, bm, eps);
            }
        } //for-loop
    }

    gThreadList.clear();

    /*if (! init_s.mPool.empty() && init_s.hasResources()) {
        std::thread thread(solveSerial, &init_s, std::ref(bm), eps);
        thread.join();
    }*/
    std::cout << "Problem is solved" << std::endl;
}

double findMin(const BM& bm, double eps, int maxstep, double b_coeff) {
    const int dim = bm.getDim();
    Box ibox;
    for (int i = 0; i < dim; i++) {
        ibox.emplace_back(bm.getBounds()[i].first, bm.getBounds()[i].second);
    }
    State s;
    s.mPool.push_back(ibox);
    s.mRecordVal = std::numeric_limits<double>::max();
    s.mMaxSteps = maxstep;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
#if 0
    solveSerial(s, bm, eps);
#else
    solve(s, bm, eps, b_coeff);
#endif
    end = std::chrono::system_clock::now();
    int mseconds = (std::chrono::duration_cast<std::chrono::microseconds> (end-start)).count();
    std::cout << "Time: " << mseconds << " microsecond\n";
    std::cout << "Time per subproblem: " << ((mseconds > 0) ? ((double) s.mSteps / (double) mseconds) : 0) << " miscroseconds." << std::endl;
    if (s.mSteps >= maxstep) {
        std::cout << "Failed to converge in " << maxstep << " steps\n";
    } else {
        std::cout << "Converged in " << s.mSteps << " steps\n";
    }
    
    std::cout << "BnB found = " << s.mRecordVal << std::endl;
    std::cout << " at x [ ";
    std::copy(s.mRecord.begin(), s.mRecord.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << "]\n";
    return s.mRecordVal;
}

bool testBench(const BM& bm) {
    constexpr double eps = 0.1;
    constexpr double b_coeff = 0.5;
    std::cout << "*************Testing benchmark**********" << std::endl;
    std::cout << bm;
    double v = findMin(bm, eps, gMaxStepsTotal, b_coeff);
    double diff = v - bm.getGlobMinY();
    if (diff > eps) {
        std::cout << "BnB failed for " << bm.getDesc() << " benchmark " << std::endl;
    }
    std::cout << "the difference is " << v - bm.getGlobMinY() << std::endl;
    std::cout << "****************************************" << std::endl << std::endl;
}

main(int argc, char** argv) {
    Benchmarks<double> tests;
    int i = 0;
    if(argc == 2) {
        i = std::atoi(argv[1]) - 1;
    }

    testBench(*tests[i]);

    /*
    for (auto bm : tests) {
        testBench(*bm);
    }*/
}
