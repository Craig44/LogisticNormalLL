
#ifndef MATRIXOPS_H_INCLUDED
#define MATRIXOPS_H_INCLUDED
#include<iostream>
using namespace std;
#include<cmath>

void getMatrix(float a[3][3])
{
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            cin>>a[i][j];
}

void transformMatrix(float matrix[3][3],float transformedMatrix[3][3])

{
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            transformedMatrix[i][j]=matrix[j][i];

}

void matrixDividebyScalar(float scalar,float inputMatrix[3][3],float outputMatrix[3][3])
{


    for(int i=0; i<3; i++)

        for(int j=0; j<3; j++)

            outputMatrix[i][j]=inputMatrix[i][j]/scalar;


}
void viewMatrix(float inv[3][3])
{
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
        {
            cout<<"\t"<<inv[i][j];
        }
        cout<<"\n";
    }
}
float cofactorwithsign(float data[3][3],int x,int y)
{
    float cofactor_v;

    cofactor_v   =
        data[(x + 1) % 3][(y + 1) % 3]
        * data[(x + 2) % 3][(y + 2) % 3]
        - data[(x + 1) % 3][(y + 2) % 3]
        * data[(x + 2) % 3][(y + 1) % 3];

    return cofactor_v;
}
float matrixDetrminant(float inputmatrix[3][3])
{
    float determinant =0;
    for (int j=0; j<3 ; j++ )
    {
        determinant+= inputmatrix[0][j]*cofactorwithsign(inputmatrix,0,j);

    }
    return determinant;


}
void multiplyMatrices(float first[3][3], float second[3][3],float output[3][3])
{
    for (int i=0; i<3 ; i++ )
        for (int j=0; j<3 ; j++ )
        {
            float element=0;
            for (int l=0; l<3 ; l++ )
                element+=first[i][l]*second[l][j];
            output[i][j]=element;
        }


}
void cofactorMatrix(float inputmatrix[3][3],float cofactormatrix [3][3])
{

    for (int i=0; i<3 ; i++ )
        for (int j=0; j<3 ; j++ )
            cofactormatrix[i][j]= cofactorwithsign(inputmatrix,i,j) ;

}


#endif // MATRIXOPS_H_INCLUDED
