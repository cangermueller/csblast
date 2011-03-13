// Copyright 2009, Andreas Biegert

#ifndef CS_PARALLEL_TEMPERING_H_
#define CS_PARALLEL_TEMPERING_H_

namespace cs {

template<class Abc, class TrainingPair>
struct ParallelTempering {
  ParallelTempering(int r, int n) : rank(r), nreplicas(n) {}

  virtual ~ParallelTempering() {}

  virtual void ReplicaExchange(double myloglike, double& mytemp) = 0;

  const int rank;        // rank of this replica
  const int nreplicas;   // total number of replicas
};

#ifdef PARALLEL

template<class Abc, class TrainingPair>
struct SingleExchParallelTempering : public ParallelTempering<Abc, TrainingPair> {

  SingleExchParallelTempering(int r, int n, double t, unsigned int s)
      : ParallelTempering<Abc, TrainingPair>(r, n),
        replicas(n),
        theta(t), ran(s) {
    for (int i = 0; i < nreplicas; ++i) replicas[i] = i;
  }

  virtual ~SingleExchParallelTempering() {}

  virtual void ReplicaExchange(double myloglike, double& mytemp) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (nreplicas < 2 || ran.doub() >= theta) return;

    const int kTagTemperatureSwap = 0;
    const int kTagLikelihoodSwap  = 1;
    const int l = ran(nreplicas - 1);  // pick random temperature level
    const int lower = replicas[l];     // rank of replica with lower temp
    const int upper = replicas[l+1];   // rank of replica with upper temp
    const double r = ran.doub();       // random float for swap decision
    double temp, loglike, alph;
    MPI_Status status;

    LOG(ERROR) << strprintf("%2d: l=%d lower=%d upper=%d", rank, l, lower, upper);

    if (rank == lower) {
      // Swap temperature and likelihood with upper replica
      MPI_Sendrecv(&mytemp, 1, MPI_DOUBLE, upper, kTagTemperatureSwap,
                   &temp , 1, MPI_DOUBLE, upper, kTagTemperatureSwap,
                   MPI_COMM_WORLD, &status);
      MPI_Sendrecv(&myloglike, 1, MPI_DOUBLE, upper, kTagLikelihoodSwap,
                   &loglike , 1, MPI_DOUBLE, upper, kTagLikelihoodSwap,
                   MPI_COMM_WORLD, &status);

      // Evaluate replicae exchange criterion
      alph = AcceptanceProb(mytemp, myloglike, temp, loglike);
      LOG(ERROR) << strprintf("%2d: LL=%10.2f T=%5.2f => LL=%10.2f T=%5.2f a=%.2f",
                              rank, myloglike, mytemp, loglike, temp, alph);
      if (r < alph) {
        LOG(ERROR) << strprintf("%2d: LL=%10.2f T=%5.2f => LL=%10.2f T=%5.2f",
                               rank, myloglike, mytemp, loglike, temp);

        mytemp = temp;   // swap temperatures
        replicas[l]   = upper;  // switch replica positions in replica table
        replicas[l+1] = lower;
      }

    } else if (rank == upper) {
      // Swap temperature and likelihood with lower replica
      MPI_Sendrecv(&mytemp, 1, MPI_DOUBLE, lower, kTagTemperatureSwap,
                   &temp, 1, MPI_DOUBLE, lower, kTagTemperatureSwap,
                   MPI_COMM_WORLD, &status);
      MPI_Sendrecv(&myloglike, 1, MPI_DOUBLE, lower, kTagLikelihoodSwap,
                   &loglike , 1, MPI_DOUBLE, lower, kTagLikelihoodSwap,
                   MPI_COMM_WORLD, &status);

      // Evaluate replicae exchange criterion
      alph = AcceptanceProb(temp, loglike, mytemp, myloglike);
      LOG(ERROR) << strprintf("%2d: LL=%10.2f T=%5.2f => LL=%10.2f T=%5.2f a=%.2f",
                              rank, myloglike, mytemp, loglike, temp, alph);
      if (r < alph) {
        LOG(ERROR) << strprintf("%2d: LL=%10.2f T=%5.2f <= LL=%10.2f T=%5.2f",
                               rank, loglike, temp, myloglike, mytemp);
        mytemp = temp;  // swap temperatures
      }
    }
    // Lower replica broadcasts new replica order to all other replicas
    MPI_Bcast(&replicas[0], nreplicas, MPI_INT, lower, MPI_COMM_WORLD);

  }

  double AcceptanceProb(double temp_i, double loklike_i,
                        double temp_j, double loklike_j) {
    double p = exp((loklike_j - loklike_i) * ((1.0 / temp_i) - (1.0 / temp_j)));
    return MIN(1.0, p);
  }

  using ParallelTempering<Abc, TrainingPair>::rank;
  using ParallelTempering<Abc, TrainingPair>::nreplicas;
  Vector<int> replicas;  // ranks of all replicas sorted by temperature
  double theta;          // probability for trying a replica exchange
  Ran ran;               // needed for random swap
};


template<class Abc, class TrainingPair>
struct SortedParallelTempering : public ParallelTempering<Abc, TrainingPair> {
  typedef std::pair<double, int> LoglikeRankPair;

  SortedParallelTempering(int r, int n)
      : ParallelTempering<Abc, TrainingPair>(r, n) {}

  virtual ~SortedParallelTempering() {}

  virtual void ReplicaExchange(double myloglike, double& mytemp) {
    const int kTagLikelihood  = 0;
    const int kTagTemperature = 1;
    Vector<double> loglikes(nreplicas);  // likelihoods in rank-order
    Vector<double> temps(nreplicas);     // temperatures in rank-order
    MPI_Status stat;

    // First, the master replica with rank 0 collects the log-likelihood and
    // temperature from all other replicas, then it broadcast these back to
    // all other replicas so that each replica has all the information it needs
    // to find out at which temperature it should run from now on.
    if (rank == 0) {
      loglikes[0] = myloglike;
      temps[0]    = mytemp;

      // Collect temperatures and likelihoods
      for (int r = 1; r < nreplicas; ++r) {
        double loglike, temp;
        MPI_Recv(&loglike, 1, MPI_DOUBLE, r, kTagLikelihood, MPI_COMM_WORLD, &stat);
        loglikes[r] = loglike;
        MPI_Recv(&temp, 1, MPI_DOUBLE, r, kTagTemperature, MPI_COMM_WORLD, &stat);
        temps[r] = temp;
      }
    } else {
      MPI_Send(&myloglike, 1, MPI_DOUBLE, 0, kTagLikelihood, MPI_COMM_WORLD);
      MPI_Send(&mytemp, 1, MPI_DOUBLE, 0, kTagTemperature, MPI_COMM_WORLD);
    }

    // Broadcast likelihoods and temperatures to all replicas
    MPI_Bcast(&loglikes[0], nreplicas, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&temps[0],    nreplicas, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    LOG(ERROR) << strprintf("%2d: loglikes unsorted = %s", rank,
                            StringifyRange(&loglikes[0],
                                           &loglikes[0] + nreplicas).c_str());
    LOG(ERROR) << strprintf("%2d: temps unsorted = %s", rank,
                            StringifyRange(&temps[0],
                                           &temps[0] + nreplicas).c_str());

    // Now each replica determines its new temperature by sorting the replicas
    // by their current likelihoods and picking the temperature such that the
    // best replicas runs at the lower temperatures.
    std::vector<LoglikeRankPair> loglike_rank_pairs;
    for (int r = 0; r < nreplicas; ++r)
      loglike_rank_pairs.push_back(std::make_pair(loglikes[r], r));
    std::sort(loglike_rank_pairs.begin(), loglike_rank_pairs.end());
    std::reverse(loglike_rank_pairs.begin(), loglike_rank_pairs.end());
    std::sort(&temps[0], &temps[0] + nreplicas);

    LOG(ERROR) << strprintf("%2d: loglikes sorted = %s", rank,
                            StringifyRange(&loglikes[0],
                                           &loglikes[0] + nreplicas).c_str());
    LOG(ERROR) << strprintf("%2d: temps sorted = %s", rank,
                            StringifyRange(&temps[0],
                                           &temps[0] + nreplicas).c_str());

    for (int i = 0; i < nreplicas; ++i) {
      if (loglike_rank_pairs[i].second == rank) {
        assert_eq(loglike_rank_pairs[i].first, myloglike);
        mytemp = temps[i];
      }
    }
  }

  using ParallelTempering<Abc, TrainingPair>::rank;
  using ParallelTempering<Abc, TrainingPair>::nreplicas;
};

#endif // PARALLEL

}  // namespace cs

#endif  // CS_PARALLEL_TEMPERING_H_
