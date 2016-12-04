#include <iostream>
#include <string.h>
#include <iomanip>
#include <armadillo>
#include <cmath>
#include <math.h>
#include <ctime>
using namespace std;


void trade(int N, int trans, arma::vec (&agents), double lambda, double alpha, double gamma){

    // Array of previous transactions
    arma::mat  cij = arma::zeros<arma::mat>(N,N);

    // Write variance to file
    //ofstream varFile;
    //varFile.open("var.dat", ios::out);

    // Arbitrary initial values to break first if-test

    double avgVarBlockOld = 1e100;
    double varVarBlockOld = 1e100;
    // Reset values after each simulation
    double cumVarBlock  = 0;
    double cumVar2Block = 0;
    int blockSize = 1e4;

    for (int i = 0; i< trans; i++){

        // Pick two random agents
        int agent_i = (int) rand() % N;
        int agent_j = (int) rand() % N;

        // Agents amount
        double m_i = agents(agent_i);
        double m_j = agents(agent_j);

        // Pick agent with favored wealth, previous transactions and i=!j
        double random_factor = (double) rand()/RAND_MAX;
        // Number of previous interactions
        int c = cij(agent_i,agent_j);
        while ( (pow(fabs(m_i-m_j),-alpha)*pow(c+1,gamma) < random_factor) || (agent_i==agent_j)){
            // Pick new agent
            agent_i = (int) rand() % N; //Alternative: Pick new both
            agent_j = (int) rand() % N;
            m_i = agents(agent_i);
            m_j = agents(agent_j);
            c = cij(agent_i,agent_j);
            random_factor = (double) rand()/RAND_MAX;

        }


        // Value of transaction between 0 and 1
        double epsilon = (double) rand()/RAND_MAX;

        // Transaction with savings
        double dm = (1-lambda)*(epsilon*m_j-(1-epsilon)*m_i);
        agents(agent_i) += dm;
        agents(agent_j) -= dm;

        // Transaction has taken place, let cij know! ij = ji
        cij(agent_i,agent_j) += 1;
        cij(agent_j,agent_i) += 1;




        // Generalized equilibrium test of variance
        double var = arma::var(agents); // Variance
        cumVarBlock += var;
        cumVar2Block += var*var;
        if (i % blockSize == 0) {
            double avgVar  = cumVarBlock  / blockSize;
            double avgVar2 = cumVar2Block / blockSize;
            // Variance of a block of variance
            double varVarBlock = avgVar2 - avgVar*avgVar;

            // Check if variance is
            if ((fabs(avgVarBlockOld - avgVar) / fabs(avgVarBlockOld)      < 0.2) &&
                (fabs(varVarBlock - varVarBlockOld) / fabs(varVarBlockOld) < 0.5)) {
                cout << "Transactions reached equilibrium at i = " << i << endl;
                break;
            } else {
                avgVarBlockOld = avgVar;
                varVarBlockOld = varVarBlock;
            }
            // Reset cumulative variance
            cumVarBlock  = 0;
            cumVar2Block = 0;

        }

        // Write variance to file
        if (i % 100 == 0) {
            //varFile << i << " " << var << endl;
        }
    }
}

void output(int N, arma::vec agents, string filename){
    ofstream ofile;
    ofile.open(filename, ios::out);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    for (int i = 0; i < N; i++){
        ofile << setw(15) << setprecision(8) << agents(i) << endl;
    }
    ofile.close();
}


int main(){
    // Change seed
    srand(time(NULL));

    string filename = "5e_0.5-2.0-2.0.dat"; // output file name
    double m0  =    100;    // Initial amount
    int N      =    1000;  // Number of agents
    int trans  =    1e7;    // Number of transactions
    int sims   =    1e3;    // Number of simulations
    double lambda = 0.5;     // Saving propensity
    double alpha  = 2.0;     // Similar wealth factor
    double gamma  = 2.0;     // Previous transactions factor
    arma::vec agents(N);    // Array of agents
    arma::vec totagents(N); // Total wealth of agents for all simulations

    // Timer variables
    clock_t start, end;
    // Start timer
    start = clock();

    // Begin simulations
    for (int i=0; i<sims;i++){

        // What simulation are we on
        cout << "Simulation " << i << endl;

        // Assign equal intial wealth
        agents.fill(m0);

        // Initiate transactions
        trade(N,trans,agents,lambda, alpha, gamma);

        // Sort and add equilibrium wealth to total
        totagents += arma::sort(agents);

    }
    // Mean value of all simulations
    agents = totagents/sims;

    // Write to file
    output(N, agents,filename);

    // Stop timer
    end = clock();

    cout << "Time elapsed: " << (double) (end - start)/CLOCKS_PER_SEC << " seconds." << endl;

    return 0;
}

