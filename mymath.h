

float mymax(float ** u , int m, int n){


	float tmp = -99999;
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            if (u[i][j] > tmp)
            	tmp = u[i][j];
        }
    }
    return tmp;

}

float mymin(float ** u , int m, int n){

	float tmp = 99999;
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            if (u[i][j] < tmp)
            	tmp = u[i][j];
        }
    }
    return tmp;

}
