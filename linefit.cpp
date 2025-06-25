#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <random>
//#include <opencv2/opencv.hpp>
#include <mutex>
std::mutex cout_mutex;

using namespace std;
//using namespace cv;

random_device rd;
mt19937 gen(rd());


pair<double, double> findMinMax(const vector<vector<double>>& vec) {
    double min = numeric_limits<double>::max();
    double max = numeric_limits<double>::min();

    for (const auto& innerVec : vec) {
        for (double value : innerVec) {
            if (value < min) {
                min = value;
            }
            if (value > max) {
                max = value;
            }
        }
    }

    return make_pair(min, max);
}

extern "C" {

// Function to convert a 2D vector of doubles to an OpenCV Mat (image)
//Mat vectorToImage(const vector<vector<double>>& data, double min, double max) {
//    int rows = data.size();
//    int cols = (rows > 0) ? data[0].size() : 0;
//
//    Mat image(rows, cols, CV_8UC1); // 8-bit single channel image
//
//    for (int i = 0; i < rows; ++i) {
//        for (int j = 0; j < cols; ++j) {
//            double val = data[i][j];
//            // Clipping the value to [0, 255] range
//            val = (val-min)/(max-min)*255;
//            image.at<uchar>(i, j) = static_cast<uchar>(val);
//        }
//    }
//
//    return image;
//}




double* linspace(double start, double stop, int num) {
    auto* result = new double[num];

    double step = (stop - start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        result[i] = start + i * step;
    }

    return result;
}


double to_start(double center, double extent){
    return center - extent/2;
}


double to_spacing(int x_size, double extent){
    return extent/x_size;
}


double linear_poly(double x, double a, double b){
    return x*a + b;
};

double inverse_linear_poly(double x, double a, double b){
    return x/a - b/a;
};

double quad_poly(double x, double a, double b, double c){
    return x*x*a + x*b + c;
};

double inverse_quad_poly(double x, double a, double b, double c){
    if (b*b-4*a*c+4*a*x > 0){
        return (-b+sqrt(b*b-4*a*c+4*a*x))/(2*a);
    } else {
        return 0;
    }
}

double cub_poly(double x, double a, double b, double c, double d){
    return x*x*x*a + x*x*b + x*c + d;
};


double inverse_cub_poly(double x, double a, double b, double c, double d)
{
    b /= a;
    c /= a;
    d = (d - x) / a;
    double q = (3.0 * c - pow(b, 2)) / 9.0;
    double r = (9.0 * b * c - 27.0 * d - 2.0 * pow(b, 3)) / 54.0;
    double disc = pow(q, 3) + pow(r, 2);
    double root_b = b / 3.0;
    vector<double> roots(3, 0);

    if (disc > 0) // One root
    {
        double s = r + sqrt(disc);
        s = ((s < 0) ? -pow(-s, (1.0 / 3.0)) : pow(s, (1.0 / 3.0)));
        double t = r - sqrt(disc);
        t = ((t < 0) ? -pow(-t, (1.0 / 3.0)) : pow(t, (1.0 / 3.0)));
        roots[0] = roots[1] = roots[2] = -root_b + s + t;
    }
    else if (disc == 0) // All roots real and at least two are equal
    {
        double r13 = ((r < 0) ? -pow(-r, (1.0 / 3.0)) : pow(r, (1.0 / 3.0)));
        roots[0] = -root_b + 2.0 * r13;
        roots[1] = roots[2] = -root_b - r13;
    }
    else // Only option left is that all roots are real and unequal
    {
        q = -q;
        double dum1 = q * q * q;
        dum1 = acos(r / sqrt(dum1));
        double r13 = 2.0 * sqrt(q);
        roots[0] = -root_b + r13 * cos(dum1 / 3.0);
        roots[1] = -root_b + r13 * cos((dum1 + 2.0 * M_PI) / 3.0);
        roots[2] = -root_b + r13 * cos((dum1 + 4.0 * M_PI) / 3.0);
    }

    sort(roots.begin(), roots.end());
//    cout << roots[1] << endl;
    return roots[1];
}


tuple<int, int, int, int> findMaxIndex(const vector<vector<vector<vector<double>>>>& vec) {
    if (vec.empty() || vec[0].empty())
        return {-1, -1, -1, -1}; // Return {-1, -1} if the vector is empty or contains empty vectors

    int max_h = 0;
    int max_i = 0;
    int max_j = 0;
    int max_k = 0;
    double max_val = vec[0][0][0][0];

    for (int h = 0; h < vec.size(); ++h) {
        for (int i = 0; i < vec[h].size(); ++i) {
            for (int j = 0; j < vec[h][i].size(); ++j) {
                for (int k = 0; k < vec[h][i][j].size(); ++k) {
                    if (vec[h][i][j][k] > max_val) {
                        max_val = vec[h][i][j][k];
                        max_h = h;
                        max_i = i;
                        max_j = j;
                        max_k = k;
                    }
                }
            }
        }
    }

    return {max_h, max_i, max_j, max_k};
}



void readCSV(const string& filename, vector<double>& column1, vector<double>& column2) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return;
    }

    string line;
    while (getline(file, line)) {
        istringstream iss(line);
        string token;
        vector<string> tokens;
        while (getline(iss, token, ',')) {
            tokens.push_back(token);
        }
        if (tokens.size() == 2) {
            try {
                double val1 = stod(tokens[0]);
                double val2 = stod(tokens[1]);
                column1.push_back(val1);
                column2.push_back(val2);
            } catch (const invalid_argument& e) {
                cerr << "Invalid argument: " << e.what() << endl;
            } catch (const out_of_range& e) {
                cerr << "Out of range: " << e.what() << endl;
            }
        } else {
            cerr << "Invalid line: " << line << endl;
        }
    }

    file.close();
}


double findMax(double* arr, int size) {
    if (size <= 0 || arr == nullptr) {
        // Handle invalid input
        return 0.0; // Return a default value or handle error appropriately
    }

    double maxVal = arr[0]; // Initialize maxVal with the first element of the array

    for (int i = 1; i < size; ++i) {
        if (arr[i] > maxVal) {
            maxVal = arr[i]; // Update maxVal if a larger value is found
        }
    }

    return maxVal;
}


pair<vector<double>, vector<double>> truncateVectors(const vector<double>& vec1, const vector<double>& vec2, double a, double b) {
    vector<double> truncatedVec1;
    vector<double> truncatedVec2;

    // Iterate through the first vector
    for (double value : vec1) {
        if (value > a && value < b) {
            truncatedVec1.push_back(value);
            truncatedVec2.push_back(value);
        }
    }

    return make_pair(truncatedVec1, truncatedVec2);
}

double interpolate_lines(double cubic_fac, double quadratic_fac, double spacing, double wl_start, vector<double> lines,
                         vector<double> compspec_x, vector<double> compspec_y){
    double j;
    double x0;
    double x1;
    double y0;
    double y1;
    double y;
    double sum = 0;
    vector<double> temp_lines(lines.size());
    for (int i = 0; i < lines.size(); ++i) {
        temp_lines[i] = inverse_cub_poly(lines[i], cubic_fac, quadratic_fac, spacing, wl_start);

        if (temp_lines[i] < compspec_x[0] ||
            temp_lines[i] > compspec_x[compspec_x.size() - 1]) { continue; }


        auto it = lower_bound(compspec_x.begin(), compspec_x.end(), temp_lines[i]);

        // 'it' points to the first element not less than temp_line
        j = distance(compspec_x.begin(), it);

        // Linear interpolation
        if (j == 0) {
            y = 0;//compspec_y[0];
        } else if (j == compspec_x.size() - 1) {
            y = 0;//compspec_y[compspec_size - 1];
        } else {
            x0 = compspec_x[j - 1];
            x1 = compspec_x[j];
            y0 = compspec_y[j - 1];
            y1 = compspec_y[j];
            y = y0 + (y1 - y0) * ((temp_lines[i] - x0) / (x1 - x0));
        }
        sum += y;
    }
    return sum;
}


double interpolate_lines_chisq(double cubic_fac, double quadratic_fac, double spacing, double wl_start, vector<double> lines,
                         vector<double> compspec_x, vector<double> compspec_y){
    double j;
    double x0;
    double x1;
    double y0;
    double y1;
    double y;
    double sum = 0;
    vector<double> temp_lines(lines.size());


    int n_sum = 0;
    for (int i = 0; i < lines.size(); ++i) {
        temp_lines[i] = inverse_cub_poly(lines[i], cubic_fac, quadratic_fac, spacing, wl_start);

        if (temp_lines[i] < compspec_x[0] ||
            temp_lines[i] > compspec_x[compspec_x.size() - 1]) { sum += 1.; continue; }


        auto it = lower_bound(compspec_x.begin(), compspec_x.end(), temp_lines[i]);

        // 'it' points to the first element not less than temp_line
        j = distance(compspec_x.begin(), it);

        // Linear interpolation
        if (j == 0) {
            y = 0;//compspec_y[0];
        } else if (j == compspec_x.size() - 1) {
            y = 0;//compspec_y[compspec_size - 1];
        } else {
            x0 = compspec_x[j - 1];
            x1 = compspec_x[j];
            y0 = compspec_y[j - 1];
            y1 = compspec_y[j];
            y = y0 + (y1 - y0) * ((temp_lines[i] - x0) / (x1 - x0));
            n_sum++;
        }
        sum += (y-1)*(y-1);
    }
    if (sum != 0) {
        return sum;
    }
    else {
        return 1000000.;
    }
}


double levyRejectionSampling(double mu, double c, normal_distribution<>& n_dist, uniform_real_distribution<>& u_dist) {
    while (true) {
        double u = u_dist(gen);
        double v = n_dist(gen);

        // Calculate candidate x
        double x_candidate = mu + c / (v * v);

        // Calculate the acceptance probability
        double p = sqrt(c / (2 * M_PI)) * exp(-c / (2 * (x_candidate - mu))) / pow(x_candidate - mu, 1.5);

        // Generate a uniform random number for acceptance decision
        double u2 = u_dist(gen);
        double u3 = u_dist(gen);

        // Accept or reject the candidate
        if (u2 <= p) {
            if (u3 > 0.5){
                return x_candidate;
            }
            else{
                return -x_candidate;
            }
        }
    }
}


void writeVectorToCSV(const std::vector<std::vector<double>>& outvec, const std::string& filename) {
    std::ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "Failed to open file " << filename << std::endl;
        return;
    }

    size_t n_samples = outvec[0].size();
    size_t n_columns = outvec.size();

    // Use a stringstream to build the CSV content in memory
    std::ostringstream oss;

    for (size_t i = 0; i < n_samples; ++i) {
        for (size_t j = 0; j < n_columns; ++j) {
            oss << std::setprecision(8) << outvec[j][i];
            if (j < n_columns - 1) {
                oss << ",";
            }
        }
        oss << "\n";
    }

    // Write the entire content to the file in one operation
    outFile << oss.str();

    outFile.close();
}

void fitlines_mkcmk(const double* compspec_x, const double* compspec_y, const double* lines,
                                                      int lines_size, int compspec_size, int n_samples, double wl_start,
                                                      double spacing, double quadratic_fac, double cubic_fac,
                                                      double wl_stepsize, double spacing_stepsize, double quad_stepsize,
                                                      double cub_stepsize, double wl_cov, double spacing_cov,
                                                      double quad_cov, double cub_cov, double acc_param, const string& outname){
    vector<double> temp_lines(lines_size);
    vector<double> compspec_x_vec(compspec_x, compspec_x + compspec_size);
    vector<double> compspec_y_vec(compspec_y, compspec_y + compspec_size);
    vector<double> lines_vec(lines, lines + lines_size);
    double this_correlation = interpolate_lines_chisq(cubic_fac, quadratic_fac, spacing, wl_start,
                                                lines_vec, compspec_x_vec, compspec_y_vec);


    const double wl_lo = wl_start-wl_cov/2;
    const double wl_hi = wl_start+wl_cov/2;
    const double spacing_lo = spacing-spacing_cov/2;
    const double spacing_hi = spacing+spacing_cov/2;
    const double quad_lo = quadratic_fac-quad_cov/2;
    const double quad_hi = quadratic_fac+quad_cov/2;
    const double cub_lo = cubic_fac-cub_cov/2;
    const double cub_hi = cubic_fac+cub_cov/2;

    double step_st;
    double step_sp;
    double step_quad;
    double step_cub;
    double step_num;
    double next_correlation;

    double nextwl;
    double nextspacing;
    double nextquad;
    double nextcub;

    double m = 1/(1-acc_param);
    double n = acc_param/(acc_param-1);
    vector<vector<double>> outvec(5, vector<double>(n_samples));

    normal_distribution<> step_dis(0, wl_stepsize);
    normal_distribution<> space_dis(0, spacing_stepsize);
    normal_distribution<> quad_dis(0, quad_stepsize);
    normal_distribution<> cub_dis(0, cub_stepsize);
    normal_distribution<> standard_normal(0, 1);
    uniform_real_distribution<> step_dist(0., 1.);
    int n_accepted = 0;

    auto start2 = chrono::high_resolution_clock::now();
    int n_burn_in = 1000000;

    for (int j = 0; j < n_samples+n_burn_in; ++j) {
        if (j % 100000 == 0 && j != 0){
            std::lock_guard<std::mutex> lock(cout_mutex);
            std::cout << "Thread " << omp_get_thread_num() 
                      << " | Step " << j 
                      << " | Acceptance rate: " 
                      << std::fixed << std::setprecision(2) 
                      << 100.0 * n_accepted / static_cast<double>(j) << "%\n";
        }

//        step_st = step_dis(gen);
//        step_sp = space_dis(gen);
//        step_quad = quad_dis(gen);
//        step_cub = cub_dis(gen);
        step_num = step_dist(gen);

        step_st =   levyRejectionSampling(0, wl_stepsize, standard_normal, step_dist);
        step_sp =   levyRejectionSampling(0, spacing_stepsize, standard_normal, step_dist);
        step_quad = levyRejectionSampling(0, quad_stepsize, standard_normal, step_dist);
        step_cub =  levyRejectionSampling(0, cub_stepsize, standard_normal, step_dist);

        nextwl = wl_start+step_st;
        nextspacing = spacing+step_sp;
        nextquad = quadratic_fac+step_quad;
        nextcub = cubic_fac+step_cub;

        if(!(wl_lo < nextwl && nextwl < wl_hi && spacing_lo < nextspacing && nextspacing < spacing_hi &&
            quad_lo < nextquad && nextquad < quad_hi && cub_lo < nextcub && nextcub < cub_hi)){
            if (j >= n_burn_in){
                outvec[0][j-n_burn_in] = wl_start;
                outvec[1][j-n_burn_in] = spacing;
                outvec[2][j-n_burn_in] = quadratic_fac;
                outvec[3][j-n_burn_in] = cubic_fac;
                outvec[4][j-n_burn_in] = this_correlation;

                //stat_outfile << setprecision(8) <<  wl_start << "," << spacing << "," << quadratic_fac << "," << cubic_fac << "," << this_correlation << "\n";
            }
            continue;
        }

        next_correlation = interpolate_lines_chisq(nextcub, nextquad, nextspacing, nextwl,
                                                   lines_vec, compspec_x_vec, compspec_y_vec);

//        cout << "Xirel " << next_correlation/this_correlation << " Triggers: " << (step_num < (20*next_correlation/this_correlation)-19) << endl;
//        cout << "This correlation: " << this_correlation << " Next correlation: " << next_correlation << endl;

        if (next_correlation < this_correlation){
            wl_start = nextwl;
            spacing = nextspacing;
            quadratic_fac = nextquad;
            cubic_fac = nextcub;
            this_correlation = next_correlation;
            n_accepted++;
        }
        else if (step_num < (m*next_correlation/this_correlation)+n){
            wl_start = nextwl;
            spacing = nextspacing;
            quadratic_fac = nextquad;
            cubic_fac = nextcub;
            this_correlation = next_correlation;
            n_accepted++;
        }

        if (j >= n_burn_in) {
            if (j == n_burn_in){
                auto end2 = chrono::high_resolution_clock::now();

                // Calculate the duration in milliseconds
                chrono::duration<double, milli> duration = end2 - start2;

            
                std::lock_guard<std::mutex> lock(cout_mutex);
                std::cout << "Thread " << omp_get_thread_num() 
                          << " | Burn-in took: " 
                          << duration.count() / 1000 << " s\n";

            }
            outvec[0][j-n_burn_in] = wl_start;
            outvec[1][j-n_burn_in] = spacing;
            outvec[2][j-n_burn_in] = quadratic_fac;
            outvec[3][j-n_burn_in] = cubic_fac;
            outvec[4][j-n_burn_in] = this_correlation;
//            stat_outfile << setprecision(8) << wl_start << "," << spacing << "," << quadratic_fac << "," << cubic_fac
//                         << "," << this_correlation << "\n";
        }
    }
    writeVectorToCSV(outvec, outname);
}


tuple <double, double, double, double> fitlines(double* compspec_x, const double* compspec_y,
                                                double* lines, int lines_size, int compspec_size,
                                                double center, double extent, double quadratic_ext,
                                                double cubic_ext, size_t c_size, size_t s_size,
                                                size_t q_size, size_t cub_size, double c_cov,
                                                double s_cov, double q_cov, double cub_cov, double zoom_fac, int n_refinement){
    cout << "Fitting lines..." << endl;
//    const size_t c_size = 100;   // 100;
//    const size_t s_size = 50;    // 50;
//    const size_t q_size = 100;   // 100;
//    const size_t cub_size = 100; // 100;

//    double c_cov = 100.;
//    double s_cov = 0.05;
//    double q_cov = 2.e-5;
//    double cub_cov = 2.5e-10;
//    cout << "fitlines reached" << endl;

    double final_c = 0;
    double final_s = 0;
    double final_q = 0;
    double final_cub = 0;

    double x0;
    double x1;
    double y0;
    double y1;
    double y;

//    cout << "Compspec_size: "<< compspec_size << endl;
//    cout << "Lines_size: "<< lines_size << endl;
//    cout << "Center: " << center<< endl;
//    cout << "Extent: " << extent<< endl;

    double dcub;
    double dq;
    double dc;
    double ds;
    vector<double> compspec_x_vec(compspec_x, compspec_x + compspec_size);

    double q_cov_frac = q_cov / q_size;
    double cub_cov_frac = cub_cov / cub_size;
    double c_cov_frac = c_cov / c_size;
    double s_cov_frac = s_cov / s_size;

    double q_cov_neghalf = -q_cov / 2;
    double cub_cov_neghalf = -cub_cov / 2;
    double c_cov_neghalf = -c_cov / 2;
    double s_cov_neghalf = -s_cov / 2;

    int j;
    vector<vector<vector<vector<double>>>> fit_vals(cub_size, vector<vector<vector<double>>>(q_size, vector<vector<double>>(c_size, vector<double>(s_size))));
    double sum = 0.;

    vector<double> temp_lines(lines_size);

//    cout << "loop reached" << endl;
    for (int n=0; n < n_refinement; n++) {
        cout << "Loop " << n+1 << " out of " << n_refinement << endl;
        #pragma omp parallel for collapse(4) schedule(dynamic) private(sum, y, y0, x0, y1, x1, dq, dcub, dc, ds, j) firstprivate(temp_lines)
        for (int cub_ind = 0; cub_ind < cub_size; ++cub_ind) {
            for (int q_ind = 0; q_ind < q_size; ++q_ind) {
                for (int c_ind = 0; c_ind < c_size; ++c_ind) {
                    for (int s_ind = 0; s_ind < s_size; ++s_ind) {
                        // cout << center+dc << " " << extent+extent*ds << endl;
                        // cout << to_start(center+dc, extent+extent*ds) << endl;
                        sum = 0.;
                        dq = linear_poly(double(q_ind), q_cov_frac, q_cov_neghalf);
                        dcub = linear_poly(double(cub_ind), cub_cov_frac, cub_cov_neghalf);
                        dc = linear_poly(double(c_ind), c_cov_frac, c_cov_neghalf);
                        ds = linear_poly(double(s_ind), s_cov_frac, s_cov_neghalf);
                        for (int i = 0; i < lines_size; ++i) {
                            //cout << extent+extent*ds << endl;
                            //cout << extent+extent*ds << endl;
                            //cout << lines[i] << ", " << cubic_ext+dcub << ", " << quadratic_ext + dq << ", " << to_spacing(compspec_size, extent + extent * ds)<< ", " << to_start(center + dc, extent + extent * ds) << endl;
                            temp_lines[i] = inverse_cub_poly(lines[i],
                                                             cubic_ext+dcub,
                                                             quadratic_ext + dq,
                                                             to_spacing(compspec_size, extent + extent * ds),
                                                             to_start(center + dc, extent + extent *
                                                                                             ds));// Get the x position from the lines array
                            if (temp_lines[i] < compspec_x[0] ||
                                temp_lines[i] > compspec_x[compspec_size - 1]) { continue; }


                            auto it = lower_bound(compspec_x_vec.begin(), compspec_x_vec.end(), temp_lines[i]);

                            // 'it' points to the first element not less than temp_line
                            j = distance(compspec_x_vec.begin(), it);

                            // Linear interpolation
                            if (j == 0) {
                                y = 0;//compspec_y[0];
                            } else if (j == compspec_size - 1) {
                                y = 0;//compspec_y[compspec_size - 1];
                            } else {
                                x0 = compspec_x[j - 1];
                                x1 = compspec_x[j];
                                y0 = compspec_y[j - 1];
                                y1 = compspec_y[j];
                                y = y0 + (y1 - y0) * ((temp_lines[i] - x0) / (x1 - x0));
                            }
//                            #pragma omp atomic
                            sum += y;
                        }
                        fit_vals[cub_ind][q_ind][c_ind][s_ind] = sum;
                    }
                }
            }
        }
        auto max_indices = findMaxIndex(fit_vals);
        int temp_ind = get<2>(max_indices);
        int temp_ind_2 = get<3>(max_indices);
        ofstream outFile("debug_"+to_string(n)+".txt");

        if (outFile.is_open()) {
            // Iterate over each row
            for (const auto& row : fit_vals) {
                // Iterate over each element in the row
                for (const auto& element : row) {
                    // Write the element to the file
                    outFile << element[temp_ind][temp_ind_2] << " ";
                }
                // Write newline character after each row
                outFile << "\n";
            }
            // Close the file
            outFile.close();
            cout << "Debug data has been written to debug_"+to_string(n)+".txt" << endl;
        } else {
            cerr << "Unable to open file!" << endl;
        }

        cout << "Indices:" << get<0>(max_indices) << "," << get<1>(max_indices) << "," << get<2>(max_indices) << "," << get<3>(max_indices) << endl;

        double d_final_c = linear_poly(double(get<2>(max_indices)), c_cov_frac, c_cov_neghalf);
        double d_final_s = linear_poly(double(get<3>(max_indices)), s_cov_frac, s_cov_neghalf);
        double d_final_q = linear_poly(double(get<1>(max_indices)), q_cov_frac, q_cov_neghalf);
        double d_final_cub = linear_poly(double(get<0>(max_indices)), cub_cov_frac, cub_cov_neghalf);

        cout << "Diffs:" << setprecision(numeric_limits<double>::digits10 + 1) << d_final_c << "," << d_final_s << "," << d_final_q << "," << d_final_cub << endl;

        final_c += d_final_c;
        final_s += d_final_s*(1+final_s);
        final_q += d_final_q;
        final_cub += d_final_cub;

        cout << "Final outputs:" << setprecision(numeric_limits<double>::digits10 + 1) << final_c << "," << final_s << "," << final_q << "," << final_cub << endl;

        center += d_final_c;
        extent *= (1 + d_final_s);
        quadratic_ext += d_final_q;
        cubic_ext += d_final_cub;

        cout << "New Centers:" << setprecision(numeric_limits<double>::digits10 + 1) << center << "," << extent << "," << quadratic_ext << "," << cubic_ext << endl;

        c_cov /= zoom_fac;
        s_cov /= zoom_fac;
        q_cov /= zoom_fac;
        cub_cov /= zoom_fac;

        q_cov_frac = q_cov / q_size;
        cub_cov_frac = cub_cov / cub_size;
        c_cov_frac = c_cov / c_size;
        s_cov_frac = s_cov / s_size;

        q_cov_neghalf = -q_cov / 2;
        cub_cov_neghalf = -cub_cov / 2;
        c_cov_neghalf = -c_cov / 2;
        s_cov_neghalf = -s_cov / 2;

    }

//    for (int i = 0; i < lines_size; ++i) {
//        cout << temp_lines[i] << ", ";
//    }
//    cout << endl;

    // Convert 2D vector to image
//    auto[min, max] = findMinMax(fit_vals[125]);
//    Mat image = vectorToImage(fit_vals[125], min, max);

    // Display the image
//    imshow("Image", image);
//    waitKey(0);
    return {final_cub, final_q, final_c, final_s};
}}


void readfile(const string& filename, double*& array, size_t& length) {
    ifstream infile(filename);
    if (!infile) {
        cerr << "Could not open the file " << filename << endl;
        array = nullptr;
        length = 0;
        return;
    }

    vector<double> values;
    string line;
    while (getline(infile, line)) {
        istringstream iss(line);
        double value;
        if (iss >> value) {
            values.push_back(value);
        } else {
            cerr << "Invalid line in file: "<< filename << " : " << line << endl;
        }
    }

    length = values.size();
    array = new double[length];
    for (size_t i = 0; i < length; ++i) {
        array[i] = values[i];
    }

    infile.close();
}

void writeOutput(const string& filename, double a, double b, double c, double d) {
    // Create an ofstream object for file output
    ofstream outFile(filename);

    // Check if the file was successfully opened
    if (!outFile) {
        cerr << "Error opening file for writing: " << filename << endl;
        return;
    }

    // Write the doubles to the file, one per line
    outFile << setprecision(numeric_limits<double>::digits10 + 1) << a << endl;
    outFile << setprecision(numeric_limits<double>::digits10 + 1) << b << endl;
    outFile << setprecision(numeric_limits<double>::digits10 + 1) << c << endl;
    outFile << setprecision(numeric_limits<double>::digits10 + 1) << d << endl;

    // Close the file
    outFile.close();

    // Optional: check if the file was successfully closed
    if (!outFile) {
        cerr << "Error closing file: " << filename << endl;
    }
}

int main(int argc, char *argv[]) {
    auto start = chrono::high_resolution_clock::now();
    if (argc != 6) {
        cerr << "ERROR: Invalid number of arguments! Please pass filenames for x,y arrays, lines and arguments.";
        return -1;
    }

    double *compspec_x = nullptr;
    double *compspec_y = nullptr;
    double *lines = nullptr;
    double *arguments = nullptr;
    size_t compspec_length = 0;
    size_t lines_length = 0;
    size_t arguments_length = 0;

    if (stoi(argv[5]) == 0) {
        readfile(argv[1], compspec_x, compspec_length);
        readfile(argv[2], compspec_y, compspec_length);
        readfile(argv[3], lines, lines_length);
        readfile(argv[4], arguments, arguments_length);

        double center = arguments[0];
        double extent = arguments[1];
        double quadratic_ext = arguments[2];
        double cubic_ext = arguments[3];

        auto [a, b, c, d] = fitlines(compspec_x, compspec_y, lines, static_cast<int>(lines_length),
                                     static_cast<int>(compspec_length),
                                     arguments[0], arguments[1], arguments[2], arguments[3],
                                     static_cast<size_t>(arguments[4]),
                                     static_cast<size_t>(arguments[5]), static_cast<size_t>(arguments[6]),
                                     static_cast<size_t>(arguments[7]),
                                     arguments[8], arguments[9], arguments[10], arguments[11], arguments[12],
                                     static_cast<int>(arguments[13]));

        a = cubic_ext + a;
        b = quadratic_ext + b;
        c = center - (extent * (1 + d)) / 2 + c;
        d = extent * (1 + d);
        writeOutput("temp/output.txt", a, b, d, c);

        // Get the ending time point
        auto end = chrono::high_resolution_clock::now();

        // Calculate the duration in milliseconds
        chrono::duration<double, milli> duration = end - start;

        cout << "Found wavelength solution in: " << duration.count() / 1000 << " s" << endl;
        return 0;
    }
    else{
        cout << "Using markov chain algorithm..." << endl;

        readfile(argv[1], compspec_x, compspec_length);
        readfile(argv[2], compspec_y, compspec_length);
        readfile(argv[3], lines, lines_length);
        readfile(argv[4], arguments, arguments_length);

        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < omp_get_num_procs(); ++i) {
            // Create a unique filename for each iteration
            string output_file = ".temp/mcmkc_output" + to_string(i) + ".txt";
            fitlines_mkcmk(compspec_x, compspec_y, lines, static_cast<int>(arguments[0]),
                           static_cast<int>(arguments[1]), static_cast<int>(arguments[2]),
                           arguments[3], arguments[4], arguments[5], arguments[6], arguments[7],
                           arguments[8], arguments[9], arguments[10], arguments[11], arguments[12], arguments[13],
                           arguments[14], arguments[15], output_file);
        }

        // Get the ending time point
        auto end = chrono::high_resolution_clock::now();

        // Calculate the duration in milliseconds
        chrono::duration<double, milli> duration = end - start;

        cout << "Found wavelength solution in: " << duration.count() / 1000 << " s" << endl;
        return 0;
    }
}
