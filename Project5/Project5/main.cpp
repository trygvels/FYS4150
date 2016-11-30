#include <iostream>
#include <string.h>
#include <iomanip>
#include <armadillo>
#include <cmath>

using namespace std;


void trade(int N, int trans, arma::vec (&agents), double lambda, double alpha){

    // Arbitrary initial values for equilibrium test
    double oldvar = 1e10;
    double oldmean1 = 1e10;

    for (int i = 0; i< trans; i++){
        // Pick agent with appoximately same
        int agent_i =
        int agent_j =
        double P = abs(m_i-m_j)^(alpha)


        // Pick two random agents
        int agent_i = (int) rand() % N;
        int agent_j = (int) rand() % N;

        // If same person picked, pick new one
        while (agent_i == agent_j) agent_j = (int) rand() % N;

        // Agents amount
        double m_i = agents(agent_i);
        double m_j = agents(agent_j);

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

        /*
        // Equilibrium check
        double mean1 = arma::mean(agents);
        double mean2 = arma::mean(dot(agents,agents));
        double var = mean2; - mean1*mean1;

        cout << "var:"<< var << endl;
        cout << i << endl;
        // If variance fluxuates less than 1 percent; Equilibrium
        if (abs(oldvar-var) < 0.0000000001*var){
            cout << "Variance does not fluxutate; Equalibrium reached" << endl;
            break;
        }

        // DOES NOT FLUCTUATE BECAUSE MONEY IS CONSERVED
        // If mean fluxuates less than 1 percent; Equilibrium
        if (abs(oldmean1-mean1) < 0.000001*mean1){
            cout << "Mean does not fluxutate; Equalibrium reached" << endl;
            break;
        }

        oldvar = var;
        oldmean1 = mean1;
        */


    }
}

void output(int N, arma::vec agents, string filename){
    ofstream ofile;
    ofile.open(filename, ios::app);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    for (int i = 0; i < N; i++){
        ofile << setw(15) << setprecision(8) << agents(i) << endl;
    }
    ofile.close();
}


int main(int argc, char *argv[]){
    string filename = "mony.dat"; // output file name
    double m0  =    1e6;    // Initial amount
    int N      =    500;  // Number of agents
    int trans  =    1e7;    // Number of transactions
    int sims   =    1;    // Number of simulations
    double lambda = 0;     // Saving propensity
    double alpha  = 0.5;
    arma::vec agents(N);    // Array of agents

    arma::vec totagents(N);
    for (int i=0; i<sims;i++){

        // What simulation are we on
        cout << "Simulation " << i << endl;

        // Assign equal intial wealth
        agents.fill(m0);

        cout << "m0:"<< m0 << endl;

        // Initiate transactions
        trade(N,trans,agents,lambda, alpha);

        // Sum all final amounts
        totagents += agents;
    }
    // Mean value of all simulations
    agents = totagents/sims;

    // Write to file
    output(N, agents,filename);




    return 0;
}

