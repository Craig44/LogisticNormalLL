#include<iostream>
using namespace std;
#include "MatrixOps.h"

int main()
{
    float inputMatrix[3][3],determinant;

    cout<<"ENTER THE ELEMENTS OF THE MATRIX:\n";
    getMatrix(inputMatrix);
    determinant=matrixDetrminant(inputMatrix);
    cout<<"THE DETERMINANT IS="<<determinant<<endl;
    if(determinant==0)
        cout<<"\nMATRIX IS NOT INVERSIBLE\n";
    else
    {
        float adjointMatrix[3][3],cofactormatrix[3][3],inversedMatrix[3][3];

        cofactorMatrix(inputMatrix,cofactormatrix);
        cout<<"THE cofactor Matrix is \n";
        viewMatrix(cofactormatrix);

        transformMatrix(cofactormatrix,adjointMatrix);
        cout<<"THE adjoint Matrix is \n";
        viewMatrix(adjointMatrix);

        matrixDividebyScalar(determinant,adjointMatrix,inversedMatrix);
        cout<<"THE output Matrix after inversion is \n";
        viewMatrix(inversedMatrix);
    }
    cout<<endl;
    float firstMatrix[3][3]
//    ={1,4,7,7,8,2,9,8,-4}
    ,secondMatrix[3][3]
//    ={2,3,-4,2,0,-5,5,1,3}
    ,outputMatrix[3][3];

    cout<< "enter first matrix"<<endl;
    getMatrix(firstMatrix);

    cout<<"enter second matrix"<<endl;
    getMatrix(secondMatrix);
    multiplyMatrices(firstMatrix,secondMatrix,outputMatrix);
    cout <<"fitst x second is "<<endl;
    viewMatrix(outputMatrix);
}
