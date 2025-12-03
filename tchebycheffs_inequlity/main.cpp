//
//  main.cpp
//  tchebycheffs_inequlity
//
//  Created by Joe on 12/2/25.
//
// supervisor: Wajeeb Gharibi



#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <curl/curl.h>

const int START_YEAR = 2015;   // or 2015
const int END_YEAR   = 2025;


// ==============================
// libcurl write callback
// ==============================

size_t WriteCallback(void* contents, size_t size, size_t nmemb, void* userp) {
    size_t totalSize = size * nmemb;
    std::string* buffer = static_cast<std::string*>(userp);
    buffer->append(static_cast<char*>(contents), totalSize);
    return totalSize;
}


// ==============================
// Fetch CSV from Alpha Vantage
// ==============================

bool fetchAlphaVantageCSV(const std::string& apiKey,
                          const std::string& symbol,
                          std::string& outCsv) {
    CURL* curl = curl_easy_init();
    if (!curl) {
        std::cerr << "Error: could not initialize CURL.\n";
        return false;
    }

    std::string url =
        "https://www.alphavantage.co/query?function=HISTORICAL_OPTIONS"
        "&symbol=" + symbol +
        "&outputsize=full"
        "&datatype=csv"
        "&apikey=" + apiKey;

    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &outCsv);

    CURLcode res = curl_easy_perform(curl);
    if (res != CURLE_OK) {
        std::cerr << "CURL error: " << curl_easy_strerror(res) << "\n";
        curl_easy_cleanup(curl);
        return false;
    }

    curl_easy_cleanup(curl);
    return true;
}


// ==============================
// Dynamic array grow helper
// ==============================

/**
 * @brief Ensures that a dynamic array has capacity for at least (needed) elements.
 *        If not, it allocates a new array with larger capacity and copies old data.
 *
 * @param arr       Reference to pointer to the dynamic array.
 * @param capacity  Reference to current capacity.
 * @param needed    Minimum required capacity.
 */
void ensureCapacity(double*& arr, int& capacity, int needed) {
    if (needed <= capacity) return;

    int newCapacity = (capacity == 0) ? 1024 : capacity * 2;
    while (newCapacity < needed) {
        newCapacity *= 2;
    }

    double* newArr = new double[newCapacity];
    for (int i = 0; i < capacity; ++i) {
        newArr[i] = arr[i];
    }

    delete[] arr;
    arr = newArr;
    capacity = newCapacity;
}


/**
 * @brief Safely extracts the year from a timestamp string.
 *
 * @param timestamp Expected format: YYYY-MM-DD
 * @return year if valid, or -1 if invalid
 */
int extractYear(const std::string& timestamp) {
    // Ensure timestamp is long enough
    if (timestamp.size() < 4)
        return -1;

    // Ensure first 4 chars are digits
    for (int i = 0; i < 4; ++i) {
        if (!std::isdigit(timestamp[i]))
            return -1;
    }

    return std::stoi(timestamp.substr(0, 4));
}


/**
 * @brief Parses Alpha Vantage CSV data and stores closing prices
 *        within a specified year range into a dynamically
 *        growing array.
 *
 * This function:
 *  - Reads CSV data line by line
 *  - Skips the header row
 *  - Extracts the timestamp and filters by the given year range
 *  - Extracts the closing price
 *  - Dynamically resizes the prices array as needed
 *
 * @param csvData     Raw CSV content returned by Alpha Vantage.
 * @param prices      Reference to a dynamic array of closing prices.
 * @param priceCount  Number of valid prices stored.
 * @param capacity    Current allocated capacity of the prices array.
 * @param startYear   Lower bound of year filter (inclusive).
 * @param endYear     Upper bound of year filter (inclusive).
 *
 * @return true if at least two prices were extracted (sufficient
 *         to compute returns), false otherwise.
 */


bool parseClosingPricesDynamic(const std::string& csvData,
                               double*& prices,
                               int& priceCount,
                               int& capacity,
                               int startYear,
                               int endYear) {

    // Reset count before parsing
    priceCount = 0;

    // Create a string stream to process the CSV line-by-line
    std::istringstream ss(csvData);
    std::string line;

    // Flag used to skip the CSV header row
    bool header = true;

    // Process each line of the CSV
    while (std::getline(ss, line)) {

        // Skip empty lines (defensive programming)
        if (line.empty())
            continue;

        // Skip the first row, which contains column headers
        if (header) {
            header = false;
            continue;
        }

        // Create a stream to parse individual columns in the row
        std::stringstream row(line);
        std::string cell;

        // --- Column 0: timestamp (YYYY-MM-DD) ---
        std::getline(row, cell, ',');
        int year = extractYear(cell);
        
        // Skip malformed or invalid timestamps
        if (year == -1)
            continue;

        // Ignore data outside the desired time window
        if (year < startYear || year > endYear)
            continue;

        // --- Skip columns: open, high, low (not needed) ---
        std::getline(row, cell, ','); // open price
        std::getline(row, cell, ','); // high price
        std::getline(row, cell, ','); // low price

        // --- Column 4: close price ---
        std::getline(row, cell, ',');
        double close = std::stod(cell);

        // Ensure the dynamic array has space for one more entry
        ensureCapacity(prices, capacity, priceCount + 1);

        // Store the closing price and increment count
        prices[priceCount++] = close;
    }

    // At least two prices are required to compute returns
    return priceCount >= 2;
}


// ==============================
// Compute returns (dynamic)
// ==============================

/**
 * @brief Computes simple daily returns from price data.
 *
 * @param prices      Array of prices.
 * @param priceCount  Number of price entries.
 * @param returns     Output pointer to returns array (will be allocated).
 * @param returnCount Output count of returns.
 */
void computeReturns(const double* prices,
                    int priceCount,
                    double*& returns,
                    int& returnCount) {
    if (priceCount < 2) {
        returns = nullptr;
        returnCount = 0;
        return;
    }

    returnCount = priceCount - 1;
    returns = new double[returnCount];

    for (int i = 1; i < priceCount; ++i) {
        double prev = prices[i - 1];
        double curr = prices[i];
        returns[i - 1] = (curr - prev) / prev;
    }
}


// ==============================
// Statistics (mean, stddev)
// ==============================

double computeMean(const double* data, int n) {
    if (n <= 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) sum += data[i];
    return sum / static_cast<double>(n);
}

double computeStdDev(const double* data, int n, double mean) {
    if (n < 2) return 0.0;
    double sumSq = 0.0;
    for (int i = 0; i < n; ++i) {
        double diff = data[i] - mean;
        sumSq += diff * diff;
    }
    double var = sumSq / static_cast<double>(n - 1); // sample variance
    return std::sqrt(var);
}


// ==============================
// Empirical tail probability
// ==============================

double empiricalTailProb(const double* data,
                         int n,
                         double mean,
                         double stddev,
                         double k) {
    if (n <= 0 || stddev == 0.0) return 0.0;

    int count = 0;
    for (int i = 0; i < n; ++i) {
        if (std::fabs(data[i] - mean) >= k * stddev) {
            ++count;
        }
    }
    return static_cast<double>(count) / static_cast<double>(n);
}


// ==============================
// main: puts everything together
// ==============================

int main(int argc, const char * argv[]) {
    // ---- user configuration ----
    std::string apiKey = "WWGYXWJZQYDH9SOT"; // <- put your key
    std::string symbol = "AAPL";
    double k = 2.0;

    const int START_YEAR = 2005;
    const int END_YEAR   = 2025;

    std::string csvData;
    std::cout << "Downloading data for symbol: " << symbol << " ...\n";

    if (!fetchAlphaVantageCSV(apiKey, symbol, csvData)) {
        std::cout << "[ERROR] fetchAlphaVantageCSV failed.\n";
        std::cout << "csvData.size() = " << csvData.size() << "\n";
        return 1;
    }

    std::cout << "Downloaded " << csvData.size() << " bytes.\n";

    // (Optional) peek at the first few hundred chars
    std::cout << "First 300 chars of response:\n"
              << csvData.substr(0, 300) << "\n\n";

    double* prices = nullptr;
    int priceCount = 0;
    int priceCapacity = 0;

    if (!parseClosingPricesDynamic(csvData,
                                   prices,
                                   priceCount,
                                   priceCapacity,
                                   START_YEAR,
                                   END_YEAR)) {
        std::cout << "[ERROR] parseClosingPricesDynamic failed.\n";
        std::cout << "priceCount = " << priceCount << "\n";
        delete[] prices;
        return 1;
    }

    std::cout << "Parsed " << priceCount << " price records.\n";

    double* returns = nullptr;
    int returnCount = 0;
    computeReturns(prices, priceCount, returns, returnCount);

    if (returnCount < 10) {
        std::cout << "[ERROR] Not enough returns to analyze. returnCount = "
                  << returnCount << "\n";
        delete[] prices;
        delete[] returns;
        return 1;
    }

    std::cout << "Sample size (returns): " << returnCount << "\n";

    double mu = computeMean(returns, returnCount);
    double sigma = computeStdDev(returns, returnCount, mu);

    std::cout << "Mean daily return:  " << mu << "\n";
    std::cout << "Std dev of returns: " << sigma << "\n";

    double chebyshevBound = 1.0 / (k * k);
    double empiricalProb  = empiricalTailProb(returns, returnCount, mu, sigma, k);

    std::cout << "\n=== Chebyshev Inequality Check ===\n";
    std::cout << "k = " << k << "\n";
    std::cout << "Chebyshev upper bound: " << chebyshevBound << "\n";
    std::cout << "Empirical probability: " << empiricalProb << "\n";

    if (empiricalProb <= chebyshevBound)
        std::cout << "Result: Inequality holds empirically.\n";
    else
        std::cout << "Result: Empirical probability exceeded the bound.\n";

    delete[] prices;
    delete[] returns;
 
    
    
    return 0;
}
