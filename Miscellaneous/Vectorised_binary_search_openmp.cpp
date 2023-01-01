#include <algorithm>
#include <iostream>
#include <vector>
#include <omp.h>
//Parallel binary search, searches for multiple values at once
int main() {
    std::vector<int> array = {1, 2, 3, 4, 5};
    std::vector<int> values = {2, 4, 6};

    // Create an array to store the results of the searches
    std::vector<int> indices(values.size(), -1);

    // Perform the binary searches in parallel using OpenMP
    #pragma omp parallel for
    for (int i = 0; i < values.size(); i++) {
        int value = values[i];
        auto it = std::lower_bound(array.begin(), array.end(), value);
        if (it != array.end() && *it == value) {
            // Value was found
            indices[i] = std::distance(array.begin(), it);
        }
    }

    // Print the results of the searches
    for (int i = 0; i < values.size(); i++)
        std::cout<<indices[i]<<" ";
}