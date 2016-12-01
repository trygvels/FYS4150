#include <iostream>
#include <string.h>
#include <iomanip>
#include <armadillo>
#include <cmath>
#include <math.h>
using namespace std;


void trade(int N, int trans, arma::vec (&agents), double lambda, double alpha, double gamma){

    // Array of previous transactions
    arma::mat  cij = arma::zeros<arma::mat>(N,N);

    // Write variance to file
    //ofstream varFile;
    //varFile.open("var.dat", ios::out);

    // Arbitrary initial values to break if-test        }

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


        /*
        // -1- Pick agent with i=!j
        while (agent_i == agent_j){
            // Pick new agent
            agent_j = (int) rand() % N;
            m_j = agents(agent_j);
        }
        */

        /*
        // -2- Pick agent with favored wealth and i=!j
        double random_factor = (double) rand()/RAND_MAX;
        while ( pow( fabs(m_i-m_j), -alpha ) > random_factor && agent_i == agent_j){
            // Pick new agent
            agent_j = (int) rand() % N;
            m_j = agents(agent_j);
        }
        */

        // -3- Pick agent with favored wealth, previous transactions and i=!j
        double random_factor = (double) rand()/RAND_MAX;
        // Number of previous interactions
        int c = cij(agent_i,agent_j)+cij(agent_j,agent_i);
        while ( (pow(fabs(m_i-m_j),-alpha)*pow(c+1,gamma) <= random_factor) || (agent_i==agent_j)){
            // Pick new agent
            agent_j = (int) rand() % N;
            m_j = agents(agent_j);
            c = cij(agent_i,agent_j)+cij(agent_j,agent_i);
            random_factor = (double) rand()/RAND_MAX;
        }
        // Value of transaction between 0 and 1
        double epsilon = (double) rand()/RAND_MAX;

        // Transaction without savings
        //double m_i_new = epsilon*(m_i+m_j);
        //double m_j_new = (1-epsilon)*(m_i+m_j);

        // Transaction with savings
        double dm = (1-lambda)*(epsilon*m_j-(1-epsilon)*m_i);
        double m_i_new = m_i + dm;
        double m_j_new = m_j - dm;

        // Amounts updated
        agents(agent_i) = m_i_new;
        agents(agent_j) = m_j_new;

        // Transaction has taken place, let cij know!
        cij(agent_i,agent_j) += 1;

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
        //    varFile << i << " " << var << endl;
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

    // Change seed for each simulation
    srand(2352276);

    string filename = "mony13.dat"; // output file name
    double m0  =    100;    // Initial amount
    int N      =    1000;  // Number of agents
    int trans  =    1e7;    // Number of transactions
    int sims   =    1e2;    // Number of simulations
    double lambda = 0.5;     // Saving propensity
    double alpha  = 1;     // Similar wealth factor
    double gamma  = 3;     // Previous transactions factor
    arma::vec agents(N);    // Array of agents
    arma::vec totagents(N); // Total wealth of agents for all simulations

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

    return 0;
}

