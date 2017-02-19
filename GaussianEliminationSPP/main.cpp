/**
 * This algorithm is used for computing gaussian elimination with
 * scaled partial pivoting. For sake of simplicity, please consider
 * vector b as part of the matrix.
 * For example, if you are dealing with a 4*4 matrix and vector b
 * which is 4*1, you simply add b to the right side of the matrix
 * and making it 4*5.
*/

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

void printMatrix(vector<vector<double>> & matrix){
    for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

vector<double> GaussianEliminationSSP(vector<vector<double>>& matrix, int n){
    //Initialize vector s
    vector<double> s;
    cout << "Initialization of s: [ ";
    for(int i = 0; i < n; i++){
        //Consider all entries except the last one, which is storing value of vector b
        double max = *max_element(matrix[i].begin(), matrix[i].end() - 1);
        cout << max << " ";
        s.push_back(max);
    }
    cout << "]" << endl;

    for(int i = 0; i < n - 1; i++){
        //Finding the row which the relative pivot element is the largest.
        double max = matrix[i][i] / s[i];
        int max_index = i;
        for(int j = i + 1; j < n; j++){
            if(max < matrix[j][i] / s[j]){
                max = matrix[j][i] / s[j];
                max_index = j;
            }
        }
        //Switching two rows, as well as two elements in vector s
        cout << "Switch row " << i+1 << " and " << max_index+1 << endl;
        swap(matrix[i], matrix[max_index]);
        swap<double>(s[i], s[max_index]);

        //Do the forward elimination
        for(int j = i + 1; j < n; j++){
            double factor = matrix[j][i] / matrix[i][i];
            matrix[j][i] = 0;
            for(int k = i + 1; k <= n; k++){
                matrix[j][k] = matrix[j][k] - factor * matrix[i][k];
            }
        }
        cout << "Forward Elimination:" << endl;
        printMatrix(matrix);

    }
    //Computing the final results
    vector<double> result;
    for (int i = n - 1; i >= 0; --i) {
        double x = matrix[i][n] / matrix[i][i];
        for(int j = 0; j < i; j++){
            matrix[j][n] -= x * matrix[j][i];
        }
        result.insert(result.begin(), x);
    }

    return result;
}

int main()
{
    int n;
    vector<vector<double>> matrix;
    cout << "Scale of the matrix:";
    cin >> n;
    cout << "Enter the entries one by one" << endl;
    for(int i = 0; i < n; i++){
        matrix.push_back({});
        for(int j = 0; j <= n; j++){
            double entry;
            cin >> entry;
            matrix[i].push_back(entry);
        }
    }

    vector<double> result = GaussianEliminationSSP(matrix, n);

    cout << "The solution of the given equations:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "x" << i+1 << "=" <<  result[i] << endl;
    }

    return 0;
}