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
#include <unordered_set>

#include <testfuncs/benchmarks.hpp>

using BM = Benchmark<double>;
using Box = std::vector<Interval<double>>;
struct State;

static int gProcs = 4;

static double gEps = 0.1;

static double gSubsSplitCoeff = 0.5;

static double gStepsSplitCoeff = 0.5;

static int gMtSubsLimit = 10;

static int gMtStepsLimit = 1000;

const static int gMaxStepsTotal = 1000000;

struct ThreadList {

   State* getReadyState();
   void addReadyState(State*);

   std::forward_list<State*>& operator()() { return mList; }
   size_t Size() const { return mThreadListSize; }
   size_t ActiveCount() const { return mActiveThreadCount; }

   void PushFront(State* st) {
       ++mThreadListSize;
       mList.push_front(st);
       //std::cout << "ThreadList size "<< mThreadListSize << std::endl;
   }

   void clear();

   std::forward_list< State* > mList;
   std::vector<State*> ready_states;

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
        mStatus(st.mStatus.load()),
        mutex() {}

    State& operator=(const State& st) {
        mRecordVal = st.mRecordVal;
        mRecord = st.mRecord;
        mPool.assign(st.mPool.begin(), st.mPool.end());
        mMaxSteps = st.mMaxSteps;
        mSteps = st.mSteps;
        mStatus.store(st.mStatus.load());

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

    bool hasResources() { return ! mPool.empty() && mMaxSteps != 0; }
    void PrintResourses() {
        std::cout << "Remaining task count: " << mPool.size() << "\n"
                     "Remainng steps count: " << mMaxSteps << "\n";
        //             "Diff: " << gMaxStepsTotal - mSteps <<std::endl;
    }

    void assignTaskTo(State* st) {
        *st = *this;
        st->mSteps = 0;
        mMaxSteps = 0;
        mPool.clear();
    }

    State* try_split(State* s = nullptr) {
        std::lock_guard<std::mutex> lock(mutex);

        const int rem_steps = mMaxSteps - mSteps;
        const int pool_size = mPool.size();

        if (pool_size <= gMtSubsLimit || rem_steps <= gMtStepsLimit) {
            return nullptr;
        }

        if (s == nullptr) {
            s = new State();
        }

        const int mv_step_count = rem_steps * gStepsSplitCoeff;

        mMaxSteps = mSteps + (rem_steps - mv_step_count);
        s->mMaxSteps = mv_step_count;

        s->mRecord = mRecord;
        s->mRecordVal = mRecordVal;

        const int mv_sub_count = pool_size * gSubsSplitCoeff;

        auto mv_beg_iter = mPool.begin() + mv_sub_count;

        s->mPool.assign(mPool.begin(), mv_beg_iter);
        mPool.erase(mPool.begin(), mv_beg_iter);

        return s;
    }

    Box getSub() {
        std::lock_guard<std::mutex> lock(mutex);
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

    //~State() { std::cout << "State destructor" << std::endl; }
    int getRemSteps() { return mMaxSteps - mSteps; }

    double mRecordVal;

    std::vector<double> mRecord;

    std::vector<Box> mPool;

    int mMaxSteps;

    int mSteps = 0;

    std::atomic<Status> mStatus;

    std::mutex mutex;
};

State* ThreadList::getReadyState() {
    if(this->ready_states.empty()) {
        return nullptr;
    }

    State* rstate = *(this->ready_states.end()-1);
    this->ready_states.pop_back();
    return rstate;
}

void ThreadList::addReadyState(State* sptr) { this->ready_states.push_back(sptr); }

void ThreadList::clear() {
    for(auto it = mList.begin(); it != mList.end();) {
        delete(*it);
        ++it;
        mList.pop_front();
    }

    this->ready_states.clear();
    mThreadListSize = 0;
    mActiveThreadCount = 0;
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

void split(const Box& ibox, State* s) {
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

    std::lock_guard<std::mutex> lock(s->mutex);
    s->mPool.push_back(std::move(b1));
    s->mPool.push_back(std::move(b2));
}

void getCenter(const Box& ibox, std::vector<double>& c) {
    const int n = ibox.size();
    for (int i = 0; i < n; i++) {
        c[i] = 0.5 * (ibox[i].lb() + ibox[i].rb());
    }
}

void solveSerial(State* s, const BM& bm) {
    const int dim = bm.getDim();
    std::vector<double> c(dim);
    while (!s->mPool.empty()) {
        Box b = s->getSub();    //Lock state's mutex
        s->mSteps ++;
        getCenter(b, c);
        double v = bm.calcFunc(c);
        if (v < s->mRecordVal) {
            s->mRecordVal = v;
            s->mRecord = c;
        }
        auto lb = bm.calcInterval(b).lb();
        if (lb <= s->mRecordVal - gEps) {
            split(b, s);    //Lock state's mutex
        }
        if (s->mSteps >= s->mMaxSteps)
            break;
    }

    s->setFinished();
}

void runThread(State* s, const BM& bm) {
    ++gThreadList.mActiveThreadCount;
    s->setProcessing();
    std::thread init_thread(solveSerial, s, std::ref(bm));
    init_thread.detach();
}

bool try_assign_run(State* sender, State* receiver, const BM& bm) {
    if (sender->hasResources() && sender->mMaxSteps > gMtStepsLimit) {
        sender->assignTaskTo(receiver);
        runThread(receiver, bm);
        return true;
    } else {
        return false;
    }
}

void solve(State& init_s, const BM& bm) {
    if (gSubsSplitCoeff >= 1) {
        std::cerr << "Warn: balance coefficient should be less than 1."
        << " gSubsSplitCoeff set to 0.5" << std::endl;
        gSubsSplitCoeff = 0.5;
    }

    State* first_s = new State();
    init_s.assignTaskTo(first_s);
    gThreadList.PushFront(first_s);

    runThread(first_s, bm);

    while(gThreadList.ActiveCount() || (init_s.hasResources() && init_s.mMaxSteps > gMtStepsLimit)) {
        for( auto iter = gThreadList().begin() ; iter != gThreadList().end(); ++iter)
        {
            State* cur_state = *iter;
            if (cur_state->isReady()) {
                continue;
            }

            if(cur_state->isFinished())
            {
                init_s.merge(*cur_state);
                cur_state->setReady();
                --gThreadList.mActiveThreadCount;
                bool flag = try_assign_run(&init_s, cur_state, bm);
                if(! flag) {
                    //std::cout << "Add ready" << std::endl;
                    gThreadList.addReadyState(cur_state);
                }
                continue;
            }

            if(gThreadList.mActiveThreadCount < gProcs)
            {
                //Load balancing procedures
                    State* new_state;
                    if (gThreadList.mActiveThreadCount == gThreadList.Size()) {
                        //std::cout << "Load balancing" << std::endl;
                        new_state = cur_state->try_split();
                        if( new_state == nullptr )
                            continue;

                        gThreadList.PushFront(new_state);
                    } else {
                        State* ready_state = gThreadList.getReadyState();
                        if( ready_state == nullptr ) {
                            std::cout << "Can't find ready. Ready size : " << gThreadList.ready_states.size() << std::endl;
                            continue;
                        }

                        new_state = cur_state->try_split(ready_state);
                        if(new_state == nullptr) {
                            gThreadList.addReadyState(ready_state);
                            continue;
                        }
                    }

                    runThread(new_state, bm);
            } else {
                std::this_thread::yield();
            }
        } //for-loop
    }

    //Processing remaining subs
    if(init_s.hasResources()) {
        solveSerial(&init_s, bm);
    }

    gThreadList.clear();

    init_s.PrintResourses();
    std::cout << "Problem is solved" << std::endl;
}

double findMin(const BM& bm) {
    const int dim = bm.getDim();
    Box ibox;
    for (int i = 0; i < dim; i++) {
        ibox.emplace_back(bm.getBounds()[i].first, bm.getBounds()[i].second);
    }
    State s;
    s.mPool.push_back(ibox);
    s.mRecordVal = std::numeric_limits<double>::max();
    s.mMaxSteps = gMaxStepsTotal;
    std::chrono::time_point<std::chrono::steady_clock> start, end;
    start = std::chrono::steady_clock::now();
#if 0
    solveSerial(&s, bm);
#else
    solve(s, bm);
#endif
    end = std::chrono::steady_clock::now();
    int mseconds = (std::chrono::duration_cast<std::chrono::microseconds> (end-start)).count();
    std::cout << "Time: " << mseconds << " microsecond\n";
    std::cout << "Time per subproblem: " << ((mseconds > 0) ? ((double) mseconds / (double) s.mSteps) : 0) << " miscroseconds." << std::endl;
    if (s.mSteps >= gMaxStepsTotal) {
        std::cout << "Failed to converge in " << gMaxStepsTotal << " steps\n";
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
    std::cout << "*************Testing benchmark**********" << std::endl;
    std::cout << "Available thread count: " << std::thread::hardware_concurrency() << std::endl;
    std::cout << "Max thread count: " << gProcs << std::endl;
    std::cout << "Max steps count: " << gMaxStepsTotal << std::endl;
    std::cout << "Eps: " << gEps << std::endl;
    std::cout << "Subs split coeffitient count: " << gSubsSplitCoeff << std::endl;
    std::cout << "Steps split coeffitient count: " << gStepsSplitCoeff << std::endl;
    std::cout << "Steps limit: " << gMtStepsLimit << std::endl;
    std::cout << "Subs limit: " << gMtSubsLimit << std::endl << std::endl;
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

main(int argc, char** argv) {
    Benchmarks<double> tests;

    std::unordered_map<std::string, int> task_map({{"Biggs EXP5 Function", 9}, {"Helical Valley function", 9},
                                            {"Hosaki function", 9}, {"Langerman-5 function", 9},
                                            {"Quintic function", 19}, {"Deckkers-Aarts function", 15},
                                            {"Dolan function", 25}, {"Goldstein Price function", 13},
                                            {"Mishra 9 function", 18}, {"Trid 10 function", 17}});

    switch (argc) {
        int temp;
        case 4: { temp = std::atoi(argv[3]); gMtStepsLimit = temp ? temp : gMtStepsLimit; }
        case 3: { temp = std::atoi(argv[2]); gProcs = temp ? temp : gProcs; }
        case 2: {int i = std::atoi(argv[1]) - 1; testBench(**(tests.begin() + i)); break;}

        default: {
            int bm_count = 0;
            int failed_tests = 0;
            std::vector<std::shared_ptr<Benchmark<double>>> failed_bm_vec;
            for (auto bm : tests) {
                std::string func_name = bm->getDesc();
                auto test_f = task_map.find(func_name);
                if(test_f == task_map.end()) {
                    continue;
                }

                gProcs = task_map[func_name];

                std::cout <<"Test " << bm_count++ << std::endl;
                int res = testBench(*bm);
                if (! res) {
                    failed_bm_vec.push_back(bm);
                    ++failed_tests;
                }
            }

            std::cout << "Failed tests count: " << failed_tests << std::endl;
            std::cout << "Failed tests:" << std::endl;
            for (auto bm : failed_bm_vec) {
                std::cout << *bm << std::endl;
            }
        }
    }
}
