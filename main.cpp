#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <fstream>

using namespace std;

int n;
vector<vector<double>> distances;
vector<vector<double>> pheromones;

const double alpha = 1.0;
const double betaa = 2.0;
const double rho = 0.05;
const double tau_max = 1.0;
const double tau_min = 0.01;
const int num_ants = 100;
const int max_iterations = 100;

struct Ant
{
    vector<int> tour;
    double tour_length;
};

vector<vector<double>> readDistanceMatrix(const string &filename, int &numCities)
{
    ifstream file(filename);
    if (!file.is_open())
    {
        cerr << "Không thể mở file: " << filename << endl;
        exit(1);
    }

    file >> numCities;
    vector<vector<double>> distMatrix(numCities, vector<double>(numCities));

    for (int i = 0; i < numCities; ++i)
    {
        for (int j = 0; j < numCities; ++j)
        {
            file >> distMatrix[i][j];
        }
    }

    file.close();
    return distMatrix;
}

void initialize_pheromones()
{
    pheromones.assign(n, vector<double>(n, tau_max));
}

double calculate_tour_length(const vector<int> &tour)
{
    double length = 0.0;
    for (int i = 0; i < n - 1; ++i)
    {
        length += distances[tour[i]][tour[i + 1]];
    }
    length += distances[tour[n - 1]][tour[0]];
    return length;
}

int select_next_city(const vector<int> &visited, int current_city)
{
    vector<double> probabilities(n, 0.0);
    double sum_probabilities = 0.0;

    for (int i = 0; i < n; ++i)
    {
        if (find(visited.begin(), visited.end(), i) == visited.end())
        {
            probabilities[i] = pow(pheromones[current_city][i], alpha) *
                               pow(1.0 / distances[current_city][i], betaa);
            sum_probabilities += probabilities[i];
        }
    }

    double r = (double)rand() / RAND_MAX;
    double cumulative_probability = 0.0;

    for (int i = 0; i < n; ++i)
    {
        if (probabilities[i] > 0)
        {
            cumulative_probability += probabilities[i] / sum_probabilities;
            if (r <= cumulative_probability)
            {
                return i;
            }
        }
    }
    return -1;
}

void update_pheromones(const vector<Ant> &ants)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            pheromones[i][j] *= (1.0 - rho);
        }
    }
    double bestTour = INT_MAX;
    for (const Ant &ant : ants)
    {
        if (bestTour > 1 / ant.tour_length)
        {
            bestTour = 1 / ant.tour_length;
        }
    }
    for (const Ant &ant : ants)
    {
        double delta_pheromone = 1.0 / ant.tour_length;
        if (delta_pheromone == bestTour)
        {
            for (int i = 0; i < n - 1; ++i)
            {
                pheromones[ant.tour[i]][ant.tour[i + 1]] += rho * tau_max;
            }
            pheromones[ant.tour[n - 1]][ant.tour[0]] += rho * tau_max;
        }
        else
        {
            for (int i = 0; i < n - 1; ++i)
            {
                pheromones[ant.tour[i]][ant.tour[i + 1]] += rho * tau_min;
            }
            pheromones[ant.tour[n - 1]][ant.tour[0]] += rho * tau_min;
        }
    }
}

vector<int> SMMAS()
{
    initialize_pheromones();

    vector<int> best_tour;
    int best_length = INT_MAX;

    for (int iter = 0; iter < max_iterations; ++iter)
    {
        vector<Ant> ants(num_ants);
        for (Ant &ant : ants)
        {
            ant.tour.push_back(rand() % n);
            while (ant.tour.size() < n)
            {
                int next_city = select_next_city(ant.tour, ant.tour.back());
                if (next_city != -1)
                {
                    ant.tour.push_back(next_city);
                }
            }
            ant.tour_length = calculate_tour_length(ant.tour);

            if (ant.tour_length < best_length)
            {
                best_length = ant.tour_length;
                best_tour = ant.tour;
            }
        }
        update_pheromones(ants);

        cout << "Iteration " << iter << ": Best tour length = " << best_length << endl;
    }

    return best_tour;
}

int main()
{
    srand(time(0));
    distances = readDistanceMatrix("n=200.txt", n);
    vector<int> best_tour = SMMAS();

    cout << "Best tour found: ";
    for (int city : best_tour)
    {
        cout << city << " ";
    }
    cout << endl;

    return 0;
}
