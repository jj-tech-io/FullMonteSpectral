//function to calculate fresnel reflection coefficient
//math
#pragma once
#include "my_functions.h"
#include <algorithm>

double getDeoxyHbValue(int wavelength)
{
    // Convert the vectors to maps for efficient access:
    std::map<int, double> deoxy_hb_map(deoxy_hb.begin(), deoxy_hb.end());
    if (deoxy_hb_map.find(wavelength) != deoxy_hb_map.end()) {
        return deoxy_hb_map[wavelength];

    }
    else {
        return 0.0;
    }
}

double getOxyHbValue(int wavelength)
{
    std::map<int, double> oxy_hb_map(oxy_hb.begin(), oxy_hb.end());

    if (oxy_hb_map.find(wavelength) != oxy_hb_map.end()) {
        return oxy_hb_map[wavelength];
    }
    else {
        return 0.0;
    }
}

double RFresnel(double n1, double n2, double cosT1) {
    double r = 0.0;
    double cosT2 = 0.0;
    const double COSZERO = 1.0 - 1.0e-12;
    const double COS90D = 1.0 * pow(10, -6);

    if (n1 == n2) { // Matched boundary
        r = 0.0;
        cosT2 = cosT1;
    }
    else if (cosT1 > COSZERO) { // Normal incident
        cosT2 = 0.0;
        r = (n2 - n1) / (n2 + n1);
        r *= r;
    }
    else if (cosT1 < COS90D) { // Very slant
        cosT2 = 0.0;
        r = 1.0;
    }
    else { // General case
        double sinT1 = std::sqrt(1 - cosT1 * cosT1);
        double sinT2 = n1 * sinT1 / n2;

        if (sinT2 >= 1.0) {
            r = 1.0;
            cosT2 = 0.0;
        }
        else {
            cosT2 = std::sqrt(1 - sinT2 * sinT2);
            double cosAP = cosT1 * cosT2 - sinT1 * sinT2;
            double cosAM = cosT1 * cosT2 + sinT1 * sinT2;
            double sinAP = sinT1 * cosT2 + cosT1 * sinT2;
            double sinAM = sinT1 * cosT2 - cosT1 * sinT2;
            r = 0.5 * sinAM * sinAM * (cosAM * cosAM + cosAP * cosAP) / (sinAP * sinAP * cosAM * cosAM);
        }
    }
    return r;
}

double xFit_1931(double wave)
{
    double t1 = (wave - 442.0) * ((wave < 442.0) ? 0.0624 : 0.0374);
    double t2 = (wave - 599.8) * ((wave < 599.8) ? 0.0264 : 0.0323);
    double t3 = (wave - 501.1) * ((wave < 501.1) ? 0.0490 : 0.0382);
    return 0.362 * exp(-0.5 * t1 * t1) + 1.056 * exp(-0.5 * t2 * t2)
        - 0.065 * exp(-0.5 * t3 * t3);
}

double yFit_1931(double wave)
{
    double t1 = (wave - 568.8) * ((wave < 568.8) ? 0.0213 : 0.0247);
    double t2 = (wave - 530.9) * ((wave < 530.9) ? 0.0613 : 0.0322);
    return 0.821 * exp(-0.5 * t1 * t1) + 0.286 * exp(-0.5 * t2 * t2);
}

double zFit_1931(double wave)
{
    double t1 = (wave - 437.0) * ((wave < 437.0) ? 0.0845 : 0.0278);
    double t2 = (wave - 459.0) * ((wave < 459.0) ? 0.0385 : 0.0725);
    return 1.217 * exp(-0.5 * t1 * t1) + 0.681 * exp(-0.5 * t2 * t2);
}

double gamma_correction(double C) {
    double abs_C = std::abs(C);
    if (abs_C > 0.0031308) {
        return 1.055 * std::pow(abs_C, 1.0 / 2.2) - 0.055;
    }
    else {
        return 12.92 * C;
    }
}
//double gamma_correction(double C) {
//    if (C > 0.0031308) {
//        return 1.055 * std::pow(C, 1.0 / 2.2) - 0.055;
//    }
//    else if (C < -0.0031308) {
//        return -1.055 * std::pow(-C, 1.0 / 2.2) + 0.055;
//    }
//    else {
//        return 12.92 * C;
//    }
//}

std::vector<double> XYZ_to_sRGB(std::vector<double> xyz, int step_size) {
    //double x = xyz[0]/step_size;
    //double y = xyz[1]/step_size;
    //double z = xyz[2]/step_size;
    double x = xyz[0] / 10;
    double y = xyz[1] / 10;
    double z = xyz[2] / 10;

    std::vector<std::vector<double>> mat3x3 = {
        {3.2406, -1.5372, -0.4986},
        {-0.9689, 1.8758, 0.0415},
        {0.0557, -0.204, 1.057}
    };

    double r = x * mat3x3[0][0] + y * mat3x3[0][1] + z * mat3x3[0][2];
    double g = x * mat3x3[1][0] + y * mat3x3[1][1] + z * mat3x3[1][2];
    double b = x * mat3x3[2][0] + y * mat3x3[2][1] + z * mat3x3[2][2];

    r = gamma_correction(r) * 255.0;
    g = gamma_correction(g) * 255.0;
    b = gamma_correction(b) * 255.0;

    // Round to 2 decimal places
    r = std::round(r * 1000.0) / 1000.0;
    g = std::round(g * 1000.0) / 1000.0;
    b = std::round(b * 1000.0) / 1000.0;

    // Vectorized version
    std::vector<double> sRGB = {r, g, b};

    return sRGB;
}

std::vector<double> Get_RGB(std::vector<double> reflectances, int step_size) {
    std::vector<double> total = { 0.0, 0.0, 0.0 };
    int index = 0;
    std::vector<double> sRGB = { 0.0, 0.0, 0.0 };
    for (double nm : wavelengths) {

        double reflectance = reflectances[index];
        double x = xFit_1931(nm) * reflectance;
        double y = yFit_1931(nm) * reflectance;
        double z = zFit_1931(nm) * reflectance;

        //XYZ to sRGB

        total[0] += x ;
        total[1] += y ;
        total[2] += z ;

        index++;
    }

    //clip values 0 to 255
    sRGB = XYZ_to_sRGB(total,step_size);

    return sRGB;
}

void WriteRowToCSV(std::ofstream& file, const std::vector<double>& row) {
    for (size_t i = 0; i < row.size(); ++i) {
        file << row[i];
        if (i < row.size() - 1) {
            file << ",";
        }
    }
    file << "\n";
}
// Assuming you have this function somewhere in your code
std::vector<double> getWavelengthList(double start, double end, double step_size, bool flag) {
    std::vector<double> array;
    for (double i = start; i <= end; i += step_size) {
        array.push_back(i);
    }
    return array;
}

void WriteHeaderToCSV(std::ofstream& file) {
    // Cm,Ch,Bm,Bh,T,sR,sG,sB
    file << "Cm,Ch,Bm,Bh,T,sR,sG,sB,";

    // Assuming step_size is known here, for example 5. Modify as needed.
    double step_size = 10;
    //std::vector<double> wavelengths = getWavelengthList(380, 780, step_size, false);
    // Convert wavelengths to strings and join with commas
    std::ostringstream oss;
    for (size_t i = 0; i < wavelengths.size(); ++i) {
        oss << wavelengths[i];

        // If this is not the last entry in the wavelengths vector, add a comma.
        if (i != wavelengths.size() - 1) {
            oss << ",";
        }
    }


    // Append to file
    file << oss.str() << std::endl; // Finish the line

    //print header
    std::cout << "Cm,Ch,Bm,Bh,T," << oss.str() << std::endl;
}

std::pair<double, double> calculate_absorption_coefficient(double wavelength) {
    // Check if wavelength matches any of the known values
    for (size_t i = 0; i < sizeof(wavelengths) / sizeof(wavelengths[0]); ++i) {
        if (wavelength == wavelengths[i]) {
            double e_HbO2 = deoxy_data[i];
            double e_Hb = oxy_data[i];
            return std::make_pair(e_HbO2, e_Hb);
        }
    }

    // If no match found, calculate using raw coefficients
    auto it = std::lower_bound(wavelengths.begin(), wavelengths.end(), wavelength);
    int i = it - wavelengths.begin();
    if (i == 0) {
        i = 1.0;
    }
    else if (i == sizeof(wavelengths) / sizeof(wavelengths[0])) {
        i = sizeof(wavelengths) / sizeof(wavelengths[0]) - 1;
    }

    // Linear interpolation for e_HbO2
    double x0 = wavelengths[i - 1], x1 = wavelengths[i];
    double y0_HbO2 = deoxy_data[i - 1], y1_HbO2 = deoxy_data[i];
    double e_HbO2 = y0_HbO2 + (y1_HbO2 - y0_HbO2) * (wavelength - x0) / (x1 - x0);

    // Linear interpolation for e_Hb
    double y0_Hb = oxy_data[i - 1], y1_Hb = oxy_data[i];
    double e_Hb = y0_Hb + (y1_Hb - y0_Hb) * (wavelength - x0) / (x1 - x0);

    return std::make_pair(e_HbO2, e_Hb);
}




