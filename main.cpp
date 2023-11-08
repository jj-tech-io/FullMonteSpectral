#include "my_functions.h"


#define Nbins 1000
#define Nbinsp1 1001
#define PI 3.1415926
#define LIGHTSPEED 2.997925e10
#define ALIVE 1
#define DEAD 0
#define THRESHOLD 0.01
#define CHANCE 0.1
#define COS90D 1e-6
#define ONE_MINUS_COSZERO 1e-12
#define COSZERO (1.0 - 1.0e-12)
#define g 0.9
#define nt 1.33


double MonteCarlo(double epi_mua, double epi_mus, double derm_mua, double derm_mus, double epidermis_thickness, int Nphotons = 1000) {

    //double ReflBin[Nbinsp1];

    double epi_albedo = epi_mus / (epi_mus + epi_mua);
    double derm_albedo = derm_mus / (derm_mus + derm_mua);
    int NR = Nbins;
    double radial_size = 2.5;
    double r = 0.0;
    // int ir = 0;
    double dr = radial_size / NR;
    //random seed
    std::vector<double> ReflBin(NR + 1, 0.0);
    srand(time(NULL));
    for (int i = 0; i < Nbinsp1; i++) {
        ReflBin[i] = 0;
    }

    for (int i_photon = 0; i_photon < Nphotons; i_photon++) {
        double W = 1.0;
        int photon_status = ALIVE;
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        double ux, uy, uz;
        double costheta, sintheta, cospsi, sinpsi, psi, uxx, uyy, uzz;
        double s, rnd;
        int it, ir;
        double mua = epi_mua;
        double mus = epi_mus;
        double albedo = epi_albedo;
        double absorb;

        // Randomly set photon trajectory to yield an isotropic source.
        costheta = 2.0 * static_cast<double>(rand()) / RAND_MAX - 1.0;
        sintheta = sqrt(1.0 - costheta * costheta);
        psi = 2.0 * PI * static_cast<double>(rand()) / RAND_MAX;
        // std::cout << "psi: " << psi << std::endl;
        ux = sintheta * cos(psi);
        uy = sintheta * sin(psi);
        uz = (fabs(costheta)); // fabs is 

        // Propagate one photon until it dies as determined by ROULETTE or reaches the surface
        it = 0;
        const int max_iterations = 10000;
        while (true) {
            it++;
            rnd = static_cast<double>(rand()) / RAND_MAX;
            // std::cout << "rnd: " << rnd << std::endl;
            while (rnd <= 0.0) {
                rnd = static_cast<double>(rand()) / RAND_MAX;

            }
            s = -log(rnd) / (mua + mus);
            x = x + (s * ux);
            y = y + (s * uy);
            z = z + (s * uz);

            if (uz < 0) {
                double s1 = fabs(z / uz);
                x = x - (s * ux);
                y = y - (s * uy);
                z = z - (s * uz);
                x = x + (s1 * ux);
                y = y + (s1 * uy);
                z = z + (s1 * uz);

                double internal_reflectance = RFresnel(1.0, nt, -uz);
                // std::cout << "internal_reflectance: " << internal_reflectance << std::endl;
                double external_reflectance = 1 - internal_reflectance;
                r = sqrt(x * x + y * y);
                ir = static_cast<int>(r / dr);
                if (ir >= NR) {
                    ir = NR;
                }
                if (ir < 0) {
                    ir = 0;
                }
                ReflBin[ir] = ReflBin[ir] + (W * external_reflectance);
                W = internal_reflectance * W;
                uz = -uz;
                x = (s - s1) * ux;
                y = (s - s1) * uy;
                z = (s - s1) * uz;
            }

            if (z < epidermis_thickness) {
                mua = epi_mua;
                mus = epi_mus;
                albedo = epi_albedo;
            }
            else {
                mua = derm_mua;
                mus = derm_mus;
                albedo = derm_albedo;
            }

            absorb = W * (1 - albedo);
            W = W - absorb;

            // Sample for costheta
            rnd = static_cast<double>(rand()) / RAND_MAX;
            if (g == 0.0) {
                costheta = 2.0 * rnd - 1.0;
            }
            else {
                double temp = (1.0 - g * g) / (1.0 - g + 2 * g * rnd);
                costheta = (1.0 + g * g - temp * temp) / (2.0 * g);
            }
            sintheta = sqrt(1.0 - costheta * costheta);

            // Sample psi
            psi = 2.0 * PI * static_cast<double>(rand()) / RAND_MAX;
            cospsi = cos(psi);
            if (psi < PI) {
                sinpsi = sqrt(1.0 - cospsi * cospsi);
            }
            else {
                sinpsi = -sqrt(1.0 - cospsi * cospsi);
            }

            if (1 - abs(uz) <= ONE_MINUS_COSZERO) {
                uxx = sintheta * cospsi;
                uyy = sintheta * sinpsi;
                uzz = costheta * copysign(uz, -1.0);
            }
            else {
                double temp = sqrt(1.0 - uz * uz);
                uxx = sintheta * (ux * uz * cospsi - uy * sinpsi) / temp + ux * costheta;
                uyy = sintheta * (uy * uz * cospsi + ux * sinpsi) / temp + uy * costheta;
                uzz = -sintheta * cospsi * temp + uz * costheta;
            }
            //update trajectory
            ux = uxx;
            uy = uyy;
            uz = uzz;
            if (W < THRESHOLD) {
                if (static_cast<double>(rand()) / RAND_MAX <= CHANCE) {
                    W = W / CHANCE;
                }
                else {
                    photon_status = DEAD;
                }
            }
            if (photon_status == DEAD || it > max_iterations) {
                break;
            }
        }
        if (i_photon >= Nphotons) {
            break;
        }
    }
    double total_reflection = 0.0;
    double max = 0.0;
    for (int i = 0; i < Nbinsp1; i++) {
        if (ReflBin[i] > max)
        {
            max = ReflBin[i]/Nphotons;
        }
        // std::cout << "ReflBin[" << each << "] = " << ReflBin[each] << std::endl;
        total_reflection += ReflBin[i] / Nphotons;

    }
    if (total_reflection > 1.0)
    {
	    	std :: cout << "total_reflection: " << total_reflection << std::endl;
            std::cout << "max: " << max << std::endl;
	}

    return total_reflection;
}
std::vector<float> generateDistribution(float minVal, float maxVal, int numSamples, double exponent = 1.0) {
    std::vector<float> values;
    float val;
    for (int i = 0; i < numSamples; ++i) {
        val = std::pow(minVal + (maxVal - minVal) * i / (numSamples - 1), exponent);
        values.push_back(round(val * 10000) / 10000); // Round to 4 decimal places
    }
    return values;
}
std::vector<double> getWavelengths(double a, double b, double s, bool print_results = false) {
    std::vector<double> result;
    for (double val = a; val <= b; val += s) {
        result.push_back(val);
    }
    if (print_results == true) {
        for (size_t i = 0; i < result.size(); i++) {
            std::cout << result[i] << std::endl;
        }
    }
    return result;
}

std::vector<double> CalculateReflectanceRow(double Cm, double Ch, double Bm, double Bh,  double T, int Nphotons = 1000) {
    // 400 to 1000nm in 10nm steps
    int step_size = 10;
    std::vector<double> reflectances(wavelengths.size());
    std::vector<double> row;
    //total
    std::vector<double> total = { 0.0, 0.0, 0.0 };
    int index = 0;
    for (double nm : wavelengths) {
        double alpha_base = 0.0244 + 8.53 * std::exp(-(nm - 154) / 66.2);
        double alpha_em = 6.6 * std::pow(10, 10) * std::pow(nm, -3.33);
        double alpha_pm = 2.9 * std::pow(10, 14) * std::pow(nm, -4.75);
        // Calculate alpha_oxy and alpha_deoxy using appropriate data source
        // double alpha_oxy = oxy_hb[index].second;
        // double alpha_deoxy = deoxy_hb[index].second;
        auto coefficients = calculate_absorption_coefficient(nm);
        double alpha_oxy = coefficients.first;
        double alpha_deoxy = coefficients.second;
        //double A = Bh*gammaOxy + (1 - Bh)*gammaDeoxy;

        int wl = static_cast<int>(nm);
        //using 
        double epidermis = Cm * (Bm * alpha_em + (1 - Bm) * alpha_pm) + (1 - Cm) * alpha_base;
        //alpha_derm_total = Ch*(Bh*alpha_oxy + (1-Bh)*alpha_deoxy) + (1-Ch)*self.get_alpha_base(wl)
        double dermis = Ch *(0.75*alpha_oxy + (0.25)*alpha_deoxy) + (1 - Ch)*alpha_base;
        //double dermis2 = Ch * (Bh * gammaOxy + (1 - Bh) * gemmaDeoxy) + (1 - Ch) * alpha_base;
        //double dermis = Ch * (A + (1 - Ch) * alpha_base);
        double scattering_epidermis = 14.74 * std::pow(nm, -0.22) + 2.22 * std::pow(10, 11) * std::pow(nm, -4.0);
        double scattering_dermis = 0.5 * scattering_epidermis;
        // Call MonteCarlo function here and store reflectance values
        double reflectance = MonteCarlo(epidermis, scattering_epidermis, dermis, scattering_dermis, T, Nphotons);

        reflectances[index] = reflectance;
        
        double x =  xFit_1931(nm);
        double y =  yFit_1931(nm);
        double z =  zFit_1931(nm);
        double XYZ [3] = {x, y, z};

        total[0] += x*reflectance;
        total[1] += y*reflectance;
        total[2] += z*reflectance;
        index++;
        //std::cout << "Total: " << total[0] << ", " << total[1] << ", " << total[2] << std::endl;
    }

    std::vector<double> sRGB = XYZ_to_sRGB(total, step_size);
    if (sRGB[0] <= 255.0 && sRGB[1] <= 255.0 && sRGB[2] <= 255.0) {
        row.push_back(Nphotons);
        row.push_back(Cm);
        row.push_back(Ch);
        row.push_back(Bm);
        row.push_back(Bh);
        row.push_back(T);
        row.insert(row.begin() + 5, sRGB.begin(), sRGB.end());
        row.insert(row.end(), reflectances.begin(), reflectances.end());
	}

    //free up memory
    total.clear();
    sRGB.clear();
    return row;
}
std::vector<double> generateSequence(double start, double end, int numSamples, double root) {
    std::vector<double> values;

    double startRootValue = std::pow(start, 1.0 / root);
    double endRootValue = std::pow(end, 1.0 / root);
    double delta = (endRootValue - startRootValue) / (numSamples - 1);

    for (int i = 0; i < numSamples; ++i) {
        double valRoot = startRootValue + delta * i;
        values.push_back(std::pow(valRoot, root));
    }

    return values;
}


std::mutex mtx; // For synchronizing output
std::mutex task_mtx; // Mutex for task queue
std::condition_variable cv; // Condition variable for the task queue

std::queue<std::function<void()>> tasks;

bool finished = false;

void ProcessAndWrite(std::ofstream& outputFile, double cm, double ch, double bm,double bh, double t, int Nphotons = 1000) {
    std::vector<double> row = CalculateReflectanceRow(cm, ch, bm,bh, t, Nphotons);
    if (row.empty()) {
            std::cout << "row is empty" << std::endl;
		    return;
	    }
    //mtx.lock();
    WriteRowToCSV(outputFile, row);
    //mtx.unlock();
    row.clear();


}
void worker() {
    while (true) {
        std::function<void()> task;

        {
            std::unique_lock<std::mutex> lock(task_mtx);

            cv.wait(lock, [] {
                return !tasks.empty() || finished;
                });

            if (tasks.empty() && finished) return;

            task = std::move(tasks.front());
            tasks.pop();
        }
        task();
    }
}
int main() {
    //clear all 
    double step_size = 10;
    int numSamples = 1;
    //std::vector<double> CmValues = generateSequence(0.00, 0.62, numSamples, 3);
    //std::vector<double> ChValues = generateSequence(0.00, 0.32, numSamples, 4);
    //std::vector<double> BmValues = generateSequence(0.00, 1.0, numSamples, 1);
    //std::vector<double> BhValues = generateSequence(0.0,1.0,numSamples,1 );
    //std::vector<double> TValues =  generateSequence(0.01,0.25, numSamples, 1);
    std::vector<double> CmValues = { 0.2 };
    std::vector<double> ChValues = { 0.1 };
    std::vector<double> BmValues = { 0.5 };
    std::vector<double> BhValues = { 0.75 };
    std::vector<double> TValues = { 0.15 };


    //   for (double val : CmValues) {
       //	std::cout << val << ", ";
       //}
    //   std::cout << std::endl;
    //   std::cout << "Ch: ";
    //   for (double val : ChValues) {
    //       std::cout << val << ", ";
    //   }
    std::cout << std::endl;

    std::cout << "size of cartesian product: " << CmValues.size() * ChValues.size() * BmValues.size() * BhValues.size() * TValues.size() << std::endl;
    //array of number of photons = [1000o to 1000000] in powers of 10
    //std::vector<int> Nphotons = { 10000, 15000, 100000};
    std::vector<int> Nphotons = { 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 110000, 120000, 130000, 140000, 150000, 160000, 170000, 180000, 190000, 200000, 210000, 220000, 230000, 240000, 250000, 260000, 270000, 280000, 290000, 300000, 310000, 320000, 330000, 340000, 350000, 360000, 370000, 380000, 390000, 400000, 410000, 420000, 430000, 440000, 450000, 460000, 470000, 480000, 490000, 500000, 510000, 520000, 530000, 540000, 550000, 560000, 570000, 580000, 590000, 600000, 610000, 620000, 630000, 640000, 650000, 660000, 670000, 680000, 690000, 700000, 710000, 720000, 730000, 740000, 750000, 760000, 770000, 780000, 790000, 800000, 810000, 820000, 830000, 840000, 850000, 860000, 870000, 880000, 890000, 900000, 910000, 920000, 930000, 940000, 950000, 960000, 970000, 980000, 990000, 1000000 };
    //array of csv file names = [converge_1000.csv, converge_10000.csv, converge_100000.csv]
    //std::vector<std::string> csvFileNames = { "converge_1000.csv", "converge_1500.csv", "converge_10000.csv", "converge_15000.csv", "converge_100000.csv" };
    std::string csvFileNames = "convergenge.csv";
    std::ofstream outputFile(csvFileNames, std::ios::out);

    //start timer
    auto start = std::chrono::high_resolution_clock::now();
    int count = 0;

    WriteHeaderToCSV(outputFile);

    const int numThreads = std::thread::hardware_concurrency();


    if (numThreads == 0) {
        std::cout << "Unable to detect number of threads. Defaulting to 4." << std::endl;
    }
    else {
        std::cout << "Detected " << numThreads << " threads." << std::endl;
    }
    std::vector<std::thread> workers;

    //for (int i = 0; i < numThreads; i++) {
    //    workers.push_back(std::thread(worker));
    //}
    for (int i = 0; i < Nphotons.size(); i++) {
        for (auto cm : CmValues) {
            for (auto ch : ChValues) {
                for (auto bm : BmValues) {
                    for (auto bh : BhValues) {
                        for (auto t : TValues) {
                            //auto task = [&, cm, ch, bm, bh, t]() {
                                ProcessAndWrite(outputFile, cm, ch, bm, bh, t, Nphotons[i]);
                                //};
                            //add a delay to allow the queue to fill up
                            /*std::this_thread::sleep_for(std::chrono::milliseconds(100));*/

                            {
                                std::unique_lock<std::mutex> lock(task_mtx);
                                //tasks.push(task);
                                //cv.notify_one();
                            }
                        }
                    }
                }
            }
        }
        std::cout << "Nphotons: " << Nphotons[i] << std::endl;
    }
    //{    
    //std::unique_lock<std::mutex> lock(task_mtx);
    //finished = true;
    //cv.notify_all(); 
    //}

        

    for (auto& worker : workers) {
        worker.join();
    }

    outputFile.close();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "elapsed time: " << elapsed.count() << " seconds" << std::endl;
    //clear CmValues, ChValues, BmValues, BhValues, TValues, Nphotons
    CmValues.clear();
    ChValues.clear();
    BmValues.clear();
    BhValues.clear();
    TValues.clear();
    Nphotons.clear();

    return 0;
}
