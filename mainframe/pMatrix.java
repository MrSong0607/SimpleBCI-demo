package mainframe;
/*
 * pMatrix.java		v0.02	May 16th, 1998
 *
 * A general purpose matrix, with 3D thingies, 
 * LU decomposition, inverse, linear system solver,
 * 
 * (05/16/1998,05/17/1998,4/2005, 5/2008, 6/2008, 
 * 7/2008,7/2009)
 *
 * Copyright(c) 1998-2009, Alex S, Particle
 */

// import java.text.*;

/*
// things to add: 

median filter
detect edges
draw line
fill polygon
texture mapping.

weighted mean
gaussian function
entropy
svd
dwt/idwt
dft/idtf (?)
sort (use 0 location as index)
percentile
median
*/


/**
 * pMatrix class, to handle all the matrix thingies.
 *
 * The 3d operations assume a 4x4 matrix.
 * Note: This thing is Right Handed!; just like OpenGL :-)
 */
public class pMatrix {

    public int rows,cols;
    
    private int pitch;
    private double[] arr;

    /**
     * update listeners.
     * these don't kick off every update. they only kick off when 
     * user calls "update"; which will update all the listeners.
     */
    private java.util.Vector listeners = null;

    public interface Listener extends java.util.EventListener {
        void update();
    }

    public void addListener(Listener o){
        if(listeners == null)
            listeners = new java.util.Vector();
        if(!listeners.contains(o))
            listeners.addElement(o);
    }
    
    public void removeListener(Listener o){
        listeners.removeElement(o);
    }

    public void removeAllListeners(){
        listeners = null;        
    }

    public void update(){
        if(listeners == null)
            return;
        for(int i=0;i<listeners.size();i++)
            ((Listener)listeners.elementAt(i)).update();
    }

    /**
     * the stack to hold working matrices (note, not very thread safe).
     * TODO: rework to make safe.
     */
    protected static java.util.Vector stack;

    /**
     * get the stack going...
     */
    static{
        stack = new java.util.Vector();
    }

    /****************************************************************
     * CONSTRUCTORS
     ****************************************************************/

    /**
     * create new matrix of size: Rows x Cols
     */
    public pMatrix(int r,int c){
        rows = r <= 0 ? 1 : r;
        cols = c <= 0 ? 1 : c;
        //pitch = (int)( Math.ceil(Math.log(cols)/Math.log(2)) );
        pitch = -1;
        while((1L<< ++pitch) < cols);
        arr = new double[rows * (1<<pitch)];
    }

    /**
     * init matrix from array
     *
     * col is whether this array needs to become a column matrix.
     */
    public pMatrix(double[] a,boolean col){
        this.rows = col ? a.length : 1;
        this.cols = col ? 1 : a.length;
        //pitch = (int)( Math.ceil(Math.log(cols)/Math.log(2)) );
        pitch = -1;
        while((1L<< ++pitch) < cols);
        arr = new double[rows * (1<<pitch)];
        System.arraycopy(a,0,arr,0,a.length);        
    }

    /**
     * make a column matrix out of an array.
     */
    public pMatrix(double[] a){
        this(a,true);
    }


    /**
     * init matrix from 2d array.
     */
    public pMatrix(double[][] a){
        this.rows = a.length;
        this.cols = a[0].length;
        //pitch = (int)( Math.ceil(Math.log(cols)/Math.log(2)) );
        pitch = -1;
        while((1L<< ++pitch) < cols);        
        arr = new double[rows * (1<<pitch)];
        for(int i=0;i<a.length;i++){
            if(a[i] != null){
                for(int j=0;j<a[i].length;j++){
                    arr[(i<<pitch)+j] = a[i][j];
                }
            }
        }
    }

    /****************************************************************
     * OTHER MATRIX CONSTRUCTORS (IDENTITY, RANDOM, ETC.)
     ****************************************************************/

    /**
     * efficient duplication 有效的复制
     */
    public pMatrix dup(){
        pMatrix v = new pMatrix(rows,cols);
        System.arraycopy(arr,0,v.arr,0,arr.length);
        return v;
    }

    /**
     * make identity matrix.
     */
    public static pMatrix I(int rows,int cols){
        pMatrix v = new pMatrix(rows,cols);
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                v.arr[(i << v.pitch) + j] = i == j ? 1.0 : 0.0;
            }
        }
        return v;
    }

    /**
     * make identity (I) matrix of same size
     */
    public pMatrix I(){
        return I(rows,cols);
    }

    /**
     * uniform random in range.均匀随机范围。
     */
    public static pMatrix urandom(int rows,int cols,double min,double max){
        pMatrix v = new pMatrix(rows,cols);
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                v.arr[(i << v.pitch) + j] = Math.random()*(max-min)+min;
            }
        }
        return v;
    }

    /**
     * uniform random from 0 to 1.
     */
    public static pMatrix urandom(int rows,int cols){
        return urandom(rows,cols,0.0,1.0);
    }

    /**
     * normal random; uses Box-Muller transform
     */
    public static pMatrix nrandom(int rows,int cols,double mean,double var){
        pMatrix v = new pMatrix(rows,cols);
        int a=0;
        double[] nr = new double[2];
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                if(a==0){   // generate?
                    // Box-Muller
                    double s = Math.sqrt(-2.0 * Math.log(Math.random()));
                    double t = 2.0*Math.PI*Math.random();
                    nr[0] = s * Math.cos(t);
                    nr[1] = s * Math.sin(t);
                }
                v.arr[(i << v.pitch) + j] = nr[a] * var + mean;
                a^=1;   // next time use next value.
            }
        }
        return v;
    }

    /**
     * normal random, mean=0, and variance (stddev) = 1
     */
    public static pMatrix nrandom(int rows,int cols){
        return nrandom(rows,cols,0.0,1.0);
    }

    /****************************************************************
     * SPECIALTY CONSTRUCTORS
     ****************************************************************/

    /**
     * returns Gram matrix; X*X^T, or Gij = xi*xj
     */
    public pMatrix gram(){
        return mult(T());
    }

    /**
     * make a Hessian matrix; Hij = yi*yj*Gij
     */
    public pMatrix hessian(pMatrix y){
        pMatrix yy = y.mult(y.T());
        return gram().smult(yy);
    }

    /****************************************************************
     * MATRIX STACK METHODS.
     ****************************************************************/

    /**
     * pushes the current matrix onto the stack,
     * (the current matrix is still there though)
     */
    public void push(){
        stack.addElement(dup());
    }

    /**
     * pops the last pushed matrix from the stack,
     * and makes it current. (the previous one is
     * erased)
     * <p>
     * NOTE: no error checking is performed, you WILL
     * get a NoSuchElementException if you're not careful
     * and try to pop an empty stack.
     *
     * @return The freshly poped matrix.
     */
    public pMatrix pop(){
        pMatrix p = (pMatrix)stack.lastElement();
        stack.removeElement(p);
        rows = p.rows;
        cols = p.cols;
        pitch = p.pitch;
        arr = p.arr;
        return this;
    }

    /****************************************************************
     * MATRIX SIZE/LAYOUT OPERATIONS.
     ****************************************************************/

    /**
     * return matrix of new size; with old data.
     */
    public pMatrix resize(int r,int c){
        pMatrix n = new pMatrix(r,c);
        r = r < rows ? r : rows;
        c = c < cols ? c : cols;
        for(int i=0;i<r;i++)
            for(int j=0;j<c;j++)
                n.set(i,j,get(i,j));
        return n;
    }

    /**
     * traverse matrix in zigzag manner and save the
     * values linearly in the output row matrix.
     *
     * assumption is that matrix is square.
     */
    public pMatrix zigzag(){
        pMatrix v = new pMatrix(rows*cols,1);
        int n = rows;
        int z = 0;
        int row = 1;
        int col = 0;
        int x,y;
        for(x = 2;x <= 2*n;x++){
            if (x <= n+1){
                y = x + 1;
                if (x % 2 == 0){ 
                    col = col + 1;
                }else{
                    row = row + 1;
                }
            }else{
                y = n+1;
                if (x%2 == 0){
                    row = row - 1; 
                    col = col + 2;
                }else{ 
                    row = row + 2; 
                    col = col - 1; 
                } 
            }
            while( (row < y) && (col < y) && (row > 0) && (col > 0)){
                v.set(z,0,get(row-1,col-1));
                z = z + 1;
                if (x%2 == 0){ 
                    row = row - 1; 
                    col = col + 1;
                }else{
                    row = row + 1; 
                    col = col - 1; 
                }
            }
        }
        return v;
    }

    /**
     * undo the zigzag operation; given an array, 
     * create a zig zag pattern.
     *
     * assumption is that matrix has enough values 
     * for a square output matrix.     
     */
    public pMatrix unzigzag(){
        int n = (int)(Math.round(Math.sqrt(rows*cols)));
        pMatrix v = new pMatrix(n,n);
        int z = 0;
        int row = 1;
        int col = 0;
        int x,y;
        for(x = 2;x <= 2*n;x++){
            if (x <= n+1){
                y = x + 1;
                if (x % 2 == 0){ 
                    col = col + 1;    
                }else{ 
                    row = row + 1; 
                }
            }else{
                y = n+1;
                if (x%2 == 0){ 
                    row = row - 1; 
                    col = col + 2;
                }else{ 
                    row = row + 2; 
                    col = col - 1; 
                }
            }
            while( (row < y) && (col < y) && (row > 0) && (col > 0)){
                v.set(row-1,col-1,get(z,0));
                z = z + 1;
                if (x%2 == 0){
                    row = row - 1;
                    col = col + 1;
                }else{ 
                    row = row + 1; 
                    col = col - 1; 
                }
            }
        }
        return v;
    }

    /****************************************************************
     * DATA ACCESS METHODS
     ****************************************************************/

    /**
     * set all to a constant
     */
    public void set(double v){
        for(int i=0;i<arr.length;i++)
            arr[i] = v;
    }

    /**
     * set at index
     */
    public void set(int i,double v){
        arr[i] = v;
    }

    /**
     * set i-row, j-column value.
     */
    public void set(int i,int j,double v){
        arr[(i << pitch) + j] = v;
    }

    /**
     * get i-row, j-column value
     */
    public double get(int i,int j){
        return arr[(i << pitch) + j];
    }

    /**
     * gets the ith value of this vector.
     * only works for 1xN or Nx1 vectors.
     */
    public double get(int i){
        return arr[i];
    }

    /**
     * get sub matrix at i-row, j-column 
     */
    public pMatrix get(int i,int j,int r,int c){
        pMatrix v = new pMatrix(r,c);
        for(int y=0;y<r;y++){
            for(int x=0;x<c;x++){
                v.set(y,x,get(i+y,j+x));
            }
        }
        return v;
    }

    /**
     * set matrix at i-row, j-column
     */
    public void set(int i,int j,pMatrix v){
        for(int y=0;y<v.rows;y++){
            for(int x=0;x<v.cols;x++){
                set(i+y,j+x,v.get(y,x));
            }
        }
    }

    /**
     * set matrix equal to passed matrix.
     */
    public void set(pMatrix v){
        set(0,0,v);
    }
    
    /**
     * turn matrix into array, left to right. 将矩阵从左到右拉成数组
     */
    public double[] toArray(){
        double[] a = new double[rows * cols];
        int k=0;
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                a[k++] = arr[(i << pitch) + j];
            }
        }
        return a;
    }

    /****************************************************************
     * BASIC MATRIX OPERATIONS (MULT, ETC.)
     ****************************************************************/

    /**
     * transverse
     */
    public pMatrix T(){
        pMatrix v = new pMatrix(cols,rows);
        if(rows == 1 || cols == 1){
            System.arraycopy(arr,0,v.arr,0,rows * cols);
        }else{
            for(int i=0;i<rows;i++){
                for(int j=0;j<cols;j++){
                    v.arr[(j << v.pitch) + i] = arr[(i << pitch) + j];
                    //v.set(j,i,get(i,j));
                }
            }
        }
        return v;
    }

    /**
     * matrix multiply.  矩阵乘法。
     * c[i][j] = sum_k a[i][k] * b[k][j]
     */
    public pMatrix mult(pMatrix b) {
        pMatrix c = new pMatrix(rows,b.cols);
        if(cols != b.rows)  // incompatible matrices.
            return null;
        for(int i=0;i<rows;i++){
            for(int j=0;j<b.cols;j++){
                double sum = 0;
                for(int k=0;k<cols;k++){
                    sum += arr[(i << pitch) + k] * b.arr[(k << b.pitch) + j];
                    //sum += get(i,k) * b.get(k,j);
                }
                c.arr[(i << c.pitch) + j] = sum;
                //c.set(i,j,sum);
            }
        }
        return c;
    }

    /**
     * matrix * vector
     */
    public double[] mult(double[] v) {
        double[] u = new double[v.length];
        for (int i=0;i<4;i++) {
            u[i] = 0;
            for (int j=0;j<4;j++)
                u[i] += arr[(i<<pitch)+j] * v[j];
        }
        return u;
    }

    /**
     * dot product between two matrices
     */
    public double dot(pMatrix b){
        return smult(b).sum();
    }

    /****************************************************************
     * SCALAR MATRIX OPERATIONS (MULT, ETC.)
     ****************************************************************/

    /**
     * scalar product; c = a * b
     */
    public pMatrix smult(pMatrix b){
        pMatrix a = this;
        pMatrix c = new pMatrix(rows,cols);
        for(int i=0;i<arr.length;i++){
            c.arr[i] = arr[i] * b.arr[i];
        }
        return c;
    }

    /**
     * scalar product; this *= b
     */
    public pMatrix smultEq(pMatrix b){
        for(int i=0;i<arr.length;i++){
            arr[i] *= b.arr[i];
        }
        return this;
    }

    /**
     * scalar product; c = a * s
     */
    public pMatrix mult(double s){
        pMatrix c = new pMatrix(rows,cols);
        for(int i=0;i<arr.length;i++){
            c.arr[i] = arr[i] * s;
        }
        return c;
    }

    /**
     * scalar product; this *= s
     */
    public pMatrix multEq(double s){
        for(int i=0;i<arr.length;i++){
            arr[i] *= s;
        }
        return this;
    }

    /**
     * scalar product; c = a / s
     */
    public pMatrix div(double s){
        pMatrix c = new pMatrix(rows,cols);
        for(int i=0;i<arr.length;i++){
            c.arr[i] = arr[i] / s;
        }
        return c;
    }

    /**
     * scalar product; this /= s
     */
    public pMatrix divEq(double s){
        for(int i=0;i<arr.length;i++){
            arr[i] /= s;
        }
        return this;
    }


    /**
     * scalar inverse, ie: c = 1/a
     */
    public pMatrix sinv(){
        pMatrix c = new pMatrix(rows,cols);
        for(int i=0;i<arr.length;i++){
            c.arr[i] = 1.0 / arr[i];
        }
        return c;        
    }

    /**
     * scalar inverse, ie: this = 1/this
     */
    public pMatrix sinvEq(){
        for(int i=0;i<arr.length;i++){
            arr[i] = 1.0 / arr[i];
        }
        return this;
    }

    /**
     * scalar sum; c = a + b
     */
    public pMatrix add(pMatrix b){
        pMatrix c = new pMatrix(rows,cols);
        for(int i=0;i<arr.length;i++){
            c.arr[i] = arr[i] + b.arr[i];
        }
        return c;
    }

    /**
     * scalar sum; this += b
     */
    public pMatrix addEq(pMatrix b){
        for(int i=0;i<arr.length;i++){
            arr[i] += b.arr[i];
        }
        return this;
    }



    /**
     * scalar sum; c = a + b
     */
    public pMatrix add(double b){
        pMatrix c = new pMatrix(rows,cols);
        for(int i=0;i<arr.length;i++){
            c.arr[i] = arr[i] + b;
        }
        return c;
    }

    /**
     * scalar sum; this += b
     */
    public pMatrix addEq(double b){
        for(int i=0;i<arr.length;i++){
            arr[i] += b;
        }
        return this;
    }


    /**
     * scalar sub; c = a - b
     */
    public pMatrix subtract(pMatrix b){
        pMatrix c = new pMatrix(rows,cols);
        for(int i=0;i<arr.length;i++){
            c.arr[i] = arr[i] - b.arr[i];
        }
        return c;
    }

    /**
     * scalar sub; this -= b
     */
    public pMatrix subtractEq(pMatrix b){
        for(int i=0;i<arr.length;i++){
            arr[i] -= b.arr[i];
        }
        return this;
    }


    /**
     * scalar sub; c = a - b
     */
    public pMatrix subtract(double b){
        pMatrix c = new pMatrix(rows,cols);
        for(int i=0;i<arr.length;i++){
            c.arr[i] = arr[i] - b;
        }
        return c;
    }

    /**
     * scalar sub; this -= b
     */
    public pMatrix subtractEq(double b){
        for(int i=0;i<arr.length;i++){
            arr[i] -= b;
        }
        return this;
    }



    /****************************************************************
     * MORE ADVANCED MATRIX OPERATIONS (SUCH AS INVERSE).
     ****************************************************************/

    /** 
     * LU decomposition
     * algorithm from: Numerical Recipes in C.
     * 
     * m is the matrix.
     * n is size of matrix (n by n).
     * pitch is (int)( Math.ceil(Math.log(n)/Math.log(2)) )
     * indx (out) is the permutation of rows (input matrix).
     *
     * the input matrix m is destroyed.
     * returns false if matrix is singular. true otherwise.
     */
    private static boolean ludcmp(double[] m,int n,int pitch,int[] indx) {
        double[] vv = new double[n];

        int i,j,k,imax=0;
        double big,temp,sum;

        for(i=0;i<n;i++){   // for every row
            big=0;
            for(j=0;j<n;j++){   // for every column
                if((temp = Math.abs(m[(i<<pitch)+j])) > big)
                    big = temp;
            }
            if(Math.abs(big) < 0.000001)    // singular matrix.
                return false;
            vv[i] = 1.0 / big;
        }

        for(j=0;j<n;j++){
            for(i=0;i<j;i++){
                sum = m[(i<<pitch)+j];
                for(k=0;k<i;k++){
                    sum -= m[(i<<pitch)+k] * m[(k<<pitch)+j];
                }
                m[(i<<pitch)+j] = sum;
            }
            big = 0;
            for(i=j;i<n;i++){
                sum = m[(i<<pitch)+j];
                for(k=0;k<j;k++)
                    sum -= m[(i<<pitch)+k] * m[(k<<pitch)+j];
                m[(i<<pitch)+j] = sum;
                if ( (temp=vv[i]*Math.abs(sum)) >= big) {
                    big=temp;
                    imax=i;
                }
            }
            if(j != imax){
                for (k=0;k<n;k++) { 
                    temp = m[(imax<<pitch)+k];
                    m[(imax<<pitch)+k] = m[(j<<pitch)+k]; 
                    m[(j<<pitch)+k] = temp;
                }
                vv[imax]=vv[j];
            }            
            indx[j]=imax;

            if (Math.abs(m[(j<<pitch)+j]) < 0.000001)
                m[(j<<pitch)+j] = 0.000001;

            if (j != n-1) { 
                temp=1.0/m[(j<<pitch)+j];
                for (i=j+1;i<n;i++)
                    m[(i<<pitch)+j] = m[(i<<pitch)+j] * temp;
            }
        }
        return true;
    }


    /** 
     * back substitution
     * algorithm from: Numerical Recipes in C.
     * uses output from ludcmp.
     * 
     * can be used to solve a system of linear equations, 
     * X a = b
     *
     * inputs are same as ludcmp, 
     * b are the input 'b' values, and output 'a' values.
     */
    private static void lubksb(double[] m,int n,int pitch,int[] indx,double[] b) {
        int i,j,ip,ii = -1;
        double sum;
        for (i=0;i<n;i++) { 
            ip=indx[i];
            sum=b[ip];
            b[ip]=b[i];
            if(ii >= 0){
                for (j=ii;j<=i-1;j++)
                    sum -= m[(i<<pitch)+j] * b[j];
            }else if(sum != 0){
                ii=i;
            }
            b[i] = sum;
        }
        for (i=n-1;i>=0;i--){ 
            sum=b[i];
            for (j=i+1;j<n;j++){ 
                sum -= m[(i<<pitch)+j] * b[j];
            }
            b[i] = sum / m[(i<<pitch)+i];
        }
    }

    /** 
     * matrix inverse
     * algorithm from: Numerical Recipes in C.
     * uses LU decomposition.
     * destroys input matrix.
     * output matrix must be the same size.
     */
    private static boolean inv(double[] m,int n,int pitch,double[] out) {

        boolean ret;
        int[] indx = new int[n];
        ret = ludcmp(m,n,pitch,indx);
        if(!ret)    // matrix is singular.
            return ret;        
        int i,j;
        double[] col = new double[n];

        for(j=0;j<n;j++){
            for(i=0;i<n;i++)
                col[i] = 0.0;
            col[j] = 1;            
            lubksb(m,n,pitch,indx,col);
            for(i=0;i<n;i++){
                out[(i<<pitch) + j] = col[i];
            }
        }
        return ret;
    }

    /**
     * get inverse of matrix.
     *
     * returns null if matrix is singular.
     */
    public pMatrix inv(){
        if(rows != cols)
            return null;
        pMatrix tmp = dup();    // we end up destroying a matrix
        pMatrix out = new pMatrix(rows,cols);
        if(inv(tmp.arr,rows,pitch,out.arr))
            return out;
        return null;
    }

    /****************************************************************
     * STATISTICAL AND AGGREGATION METHODS
     ****************************************************************/

    /**
     * trace operation; sum of diagonal
     */
    public double trace(){
        double sum = 0;
        for(int i=0;i<rows && i<cols;i++)
            sum += get(i,i);
        return 0;
    }

    /**
     * Frobenius norm, M.N = trace(M^T*N)
     */
    public double frobenius(pMatrix b){
        return T().mult(b).trace();
    }

    /**
     * similarity
     */
    public double alignment(pMatrix b){
        return frobenius(b) / Math.sqrt( frobenius(this) * b.frobenius(b) );
    }

    /**
     * distance squared
     */
    public double distSq(pMatrix b){
        double d = 0;
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                double v = get(i,j) - b.get(i,j);
                d += v*v;
            }
        }
        return d;
    }

    /**
     * euclidean distance   欧几里得距离
     */
    public double dist(pMatrix b){
        return Math.sqrt(distSq(b));
    }

    /**
     * correlations
     * returns array with: 
     *  pearson, leeandlee, chandra, coxluminiandmaio, manjunath
     */
    public static double[] correlations(pMatrix x,pMatrix y){
        double sumx,sumy,sumxy,sumx2,sumy2;
        sumx = sumy = sumxy = sumx2 = sumy2 = 0.0;
        for(int i=0;i<x.rows;i++){
            for(int j=0;j<x.cols;j++){
                double xv = x.get(i,j);
                double yv = y.get(i,j);                
                sumx  += xv;
                sumy  += yv;
                sumxy += xv*yv;
                sumx2 += xv*xv;    
                sumy2 += yv*yv;
            }
        }
        double n = x.rows * x.cols;
        double pearson = ( (n *sumxy) - (sumx * sumy) ) / 
            Math.sqrt( (n*sumx2 - sumx*sumx) * (n*sumy2 - sumy*sumy) );
        double leeandlee = sumxy / sumx2;
        double chandra = sumxy / (Math.sqrt(sumx2) * Math.sqrt(sumy2));
        double coxluminiandmaio = sumxy / Math.sqrt(sumy2);
        double manjunath = sumxy / sumy2;
        return new double[] { pearson, leeandlee, chandra, coxluminiandmaio, manjunath };
    } 

    /**
     * returns Pearson correlation
     */
    public double pearson(pMatrix b){
        double[] a = correlations(this,b);
        return a[0];
    }

    /**
     * calculate stats for this matrix.  计算该矩阵的属性
     */
    public double[] stats(){
        double sum = 0, sum2=0;
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                double v = arr[(i << pitch) + j];
                sum += v;
                sum2 += v*v;
            }
        }
        double n = rows*cols;
        double avg = sum / n;
        double variance = sum2 - (sum*sum/n);
        double stddev = Math.sqrt(variance);
        return new double[] { sum, avg, sum2, variance, stddev };  //总和，平局值，平方和，方差，标准差
    }
    
    public double ValueOfRows(int row,int n,int m){
    	row=row-1;
    	double sum =0;      
            for(int j=n-1;j<m-1;j++){
                double v = arr[(row << pitch) + j];
                sum += v; 
            }
        
    	return sum;
    }
    

    /**
     * sum up all the values in the matrix. 求和
     */
    public double sum(){
        double[] s = stats();
        return s[0];
    }

    /**
     * avg all values in matrix.  矩阵平均值
     */
    public double avg(){
        double[] s = stats();
        return s[1];
    }

    /**
     * avg all values in matrix.
     */
    public double mean(){
        return avg();
    }

    /**
     * variance
     * second moment - first moment squared/n   方差
     */
    public double variance(){
        double[] s = stats();
        return s[3];
    }

    /**
     * standard deviation  标准偏差
     */
    public double stddev(){
        double[] s = stats();
        return s[4];
    }

    /**
     * find min/max values.  找到最大值或最小值
     */
    public double[] minmax(){
        double min,max;
        min=arr[0]; max=arr[0];
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                double a = arr[(i<<pitch)+j];
                if(a < min)
                    min = a;
                if(a > max)
                    max = a;
            }
        }
        return new double[]{ min, max};
    }

    public double min(){
        double[] m = minmax();
        return m[0];
    }
    public double max(){
        double[] m = minmax();
        return m[1];
    }
    

    /**
     * magnitude (length)  大小
     */
    public double mag(){
        return Math.sqrt(dot(this));
    }

    /****************************************************************
     * NORMALIZATION METHODS
     ****************************************************************/

    /**
     * normalize the whole thing.
     */
    public pMatrix norm(){
        return mult(1/mag());
    }

    /**
     * normalize the whole thing.
     */
    public pMatrix normEq(){
        return multEq(1/mag());
    }

    /**
     * normalize rows  规范化
     */
    public pMatrix normRows(){
        pMatrix c = new pMatrix(rows,cols);    
        for(int i=0;i<rows;i++){
            double l = 0;
            for(int j=0;j<cols;j++){
                double v = get(i,j);
                l += v*v;
            }
            l = Math.sqrt(l);
            for(int j=0;j<cols;j++){
                c.set(i,j, get(i,j)/l);
            }
        }
        return c;
    }

    /**
     * normalize rows in place.
     */
    public pMatrix normRowsEq(){
        for(int i=0;i<rows;i++){
            double l = 0;
            for(int j=0;j<cols;j++){
                double v = get(i,j);
                l += v*v;
            }
            l = Math.sqrt(l);
            for(int j=0;j<cols;j++){
                set(i,j, get(i,j)/l);
            }
        }
        return this;
    }


    /****************************************************************
     * MATRIX TRANSFORMS (DCT, ETC.)
     ****************************************************************/

    /**
     * returns constants needed for DCT/IDCT functions.
     * TODO: figure out some way to cache these in statics.
     */
    private pMatrix dctC(){
        pMatrix c = new pMatrix(rows,cols);
        for(int j=0;j<cols;j++){
            c.arr[j] = 1.0/Math.sqrt(cols);
        }
        for(int i=1;i<rows;i++){
            for(int j=0;j<cols;j++){
                c.arr[(i << c.pitch) + j] = 
                    Math.sqrt(2.0 / cols) * 
                    Math.cos( Math.PI * ( 2.0 * j + 1) * i / ( 2.0 * rows));
            }
        }
        return c;
    }

    /**
     * forward DCT (discrete cosine transform).
     */
    public pMatrix dct(){
        pMatrix c = dctC();
        pMatrix ct = c.T();
        return c.mult(mult(ct));        
    }
    
    /**
     * inverse DCT (discrete cosine transform)
     */
    public pMatrix idct(){
        pMatrix c = dctC();
        pMatrix ct = c.T();
        return ct.mult(mult(c));
    }


    /****************************************************************
     * THRESHOLD FUNCTIONS, USEFUL FOR NEURAL NETWORKS.
     ****************************************************************/

    private double sigmoidfunc(double x,double c){
        return 1.0/(1.0+Math.exp(-x*c));
    }
    // first derivative of sigmoid.
    private double sigmoiddxfunc(double x,double c){
        return sigmoidfunc(x,c)*(1 - sigmoidfunc(x,c));
    }


    /**
     * Neural Network thing.
     * for every element, apply sigmoid function.
     * f(a) = 1 / (1 + exp(-a)) 
     */
    public pMatrix sigmoid(double s){
        pMatrix out = new pMatrix(rows,cols);
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                out.set(i,j, sigmoidfunc(get(i,j),s));
            }
        }
        return out;
    }

    /**
     * Neural Network thing.
     * for every element, apply sigmoid function.
     * f(a) = sigmoid(x) * (1-sigmoid(x)) 
     */
    public pMatrix sigmoiddx(double s){
        pMatrix out = new pMatrix(rows,cols);
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                out.set(i,j, sigmoiddxfunc(get(i,j),s));
            }
        }
        return out;
    }

    /**
     * Neural Network thing.
     * for every element, apply sigmoid function.
     * f(a) = 1 / (1 + exp(-a)) 
     */
    public pMatrix sigmoidEq(double s){
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                set(i,j,sigmoidfunc(get(i,j),s));
            }
        }
        return this;
    }

    /**
     * Neural Network thing.
     * for every element, apply sigmoid function.
     * f(a) = 
     */
    public pMatrix sigmoiddxEq(double s){
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                set(i,j,sigmoiddxfunc(get(i,j),s));
            }
        }
        return this;
    }

    /**
     * Neural Network thing.
     * for every element, apply sigmoid function.
     * f(a) = 1 / (1 + exp(-a)) 
     */
    public pMatrix sigmoid(){
        return sigmoid(1.0);
    }

    /**
     * Neural Network thing.
     * for every element, apply sigmoid function.
     */
    public pMatrix sigmoiddx(){
        return sigmoiddx(1.0);
    }


    /**
     * Neural Network thing.
     * for every element, apply sigmoid function.
     * f(a) = 1 / (1 + exp(-a)) 
     */
    public pMatrix sigmoidEq(){
        return sigmoidEq(1.0);
    }

    /**
     * Neural Network thing.
     * for every element, apply sigmoid function.
     * f(a) = 1 / (1 + exp(-a)) 
     */
    public pMatrix sigmoiddxEq(){
        return sigmoidEq(1.0);
    }

    /**
     * Neural Network thing
     * for every element, apply tanh (sigmoid with -1 to 1 range).
     * f(a) = tanh(a)
     */
    public pMatrix tanh(){
        pMatrix out = new pMatrix(rows,cols);
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                //double a = get(i,j);
                //double e1 = Math.exp(a);
                //double e2 = Math.exp(-a);
                //out.set(i,j,(e1-e2)*(e1+e2));
                out.set(i,j, 2.0/(1.0+Math.exp(-get(i,j)))-1);
            }
        }
        return out;
    }

    /**
     * Neural Network thing
     * for every element, apply threshold with probability sigmoid.
     * f(a) = 1 with probability 1 / (1 + exp(-a)), else -1.
     */
    public pMatrix heatbath(){
        pMatrix out = new pMatrix(rows,cols);
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                double a = 1.0/(1.0+Math.exp(-get(i,j)));
                out.set(i,j,Math.random() < a ? 1 : -1);
            }
        }
        return out;
    }
    
    /**
     * Neural Network thing.
     * apply threshold funtion
     * f(a) = a gt 0 ? 1 : -1
     */
    public pMatrix threshold(){
        pMatrix out = new pMatrix(rows,cols);
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                double a = get(i,j);
                out.set(i,j, a > 0.0 ? 1.0 : -1.0);
            }
        }
        return out;
    }

    /****************************************************************
     * VITERBI ALGO: TODO
     ****************************************************************/

    /**
     * returns probabilty of observations,
     * lp is the most likely state sequence (integers index into states)
     */
     /*
    public static double viterbi(pMatrix sp, pMatrix tp,pMatrix ep,int[] lp){

*/
/*
sub forward_viterbi {
   my ($y, $X, $sp, $tp, $ep) = @_;
   my $T = {};
   for my $state (@$X){
       ##          prob.      V. path  V. prob.
       $T->{$state} = [$sp->{$state}, [$state], $sp->{$state}];
   }
   for my $output (@$y){
       my $U = {};
       for my $next_state (@$X){
           my $total = 0;
           my $argmax;
           my $valmax = 0;
           for my $source_state (@$X){
               my ($prob, $v_path, $v_prob) = @{ $T->{$source_state} };
               my $p = $ep->{$source_state}{$output} * 
                    $tp->{$source_state}{$next_state};
               $prob *= $p;
               $v_prob *= $p;
               $total += $prob;
               if($v_prob > $valmax){
                   $argmax = [@$v_path,$next_state];
                   $valmax = $v_prob;
               }
           }
           $U->{$next_state} = [$total, $argmax, $valmax];
       }
       $T = $U;
   }
   ## apply sum/max to the final states:
   my $total = 0;
   my $argmax;
   my $valmax = 0;
   for my $state (@$X){
       my ($prob, $v_path, $v_prob) = @{$T->{$state}};
       $total += $prob;
       if($v_prob > $valmax){
           $argmax = $v_path;
           $valmax = $v_prob;
       }
   }
   return ($total, $argmax, $valmax);
}
*/

    /****************************************************************
     * LINEAR SYSTEM SOLUTIONS.
     ****************************************************************/

    /**
     * solve X * w = y
     * 
     * returns false if matrix is singular
     * w is a row matrix.
     */
    public boolean solv(double[] y){
        if(rows <= cols || cols != y.length)
            return false;
        int n = cols;
        // pick top n * n matrix.
        pMatrix tmp = get(0,0,n,n);        
        boolean ret;
        int[] indx = new int[n];
        ret = ludcmp(tmp.arr,n,pitch,indx);
        if(!ret)    // matrix is singular.
            return ret;
        lubksb(tmp.arr,n,pitch,indx,y);
        return true;
    }


    /**
     * solve X * w = y
     * 
     * y and w are column matrices.
     * returns 'w' matrix.
     */
    public pMatrix solv(pMatrix y){
        if(rows <= cols || cols != y.rows || y.cols != 1)
            return null;
        int n = cols;
        // pick top n * n matrix.
        pMatrix tmp = get(0,0,n,n);
        boolean ret;
        int[] indx = new int[n];
        ret = ludcmp(tmp.arr,n,pitch,indx);
        if(!ret)    // matrix is singular.
            return null;
        pMatrix out = y.dup();
        lubksb(tmp.arr,n,pitch,indx,out.arr);
        return out;
    }

    /****************************************************************
     * FITTING POLYNOMIALS.
     ****************************************************************/

    /**
     * fit a polynomial to a buncha data points.
     * xy is a matrix of [x,y] rows.
     * n is the degree of polynomial (number of elements in returned array).
     *
     * returns an n degree polynomial. line is a 2 degree polynomial.
     */
    public static double[] fitPoly(pMatrix xy,int n){

        // make the X (samples) matrix.
        pMatrix X = new pMatrix(xy.rows,n+1);
        pMatrix Y = new pMatrix(xy.rows,1);
        for(int i=0;i<xy.rows;i++){
            double x = xy.get(i,0);
            for(int j=0;j<=n;j++){
                X.set(i,j,Math.pow(x,j));
            }
            Y.set(i,0,xy.get(i,1));
        }

        // Ridge Regression
        //
        // below are two identical ways of solving this: primal and dual.
        // the primal method is faster when dimensions are small
        // dual method is faster when number of dimensions is high.
        // dual also doesn't require actual inputs; but an inner (dot) product
        // between input points (allowing the use of kernel functions).
        
        // primal solution
        // w = (X^T * X + aI)^-1 * X^T * y
        pMatrix Xt = X.T();
        pMatrix G = Xt.mult(X);
        pMatrix w = G.add( G.I().mult(0.01) ).inv().mult(Xt).mult(Y);
                
        /*
        // dual solution
        // alpha = (X * X^T + aI)^-1 * y
        // w = X^T * alpha
        pMatrix Xt = X.T();
        pMatrix G = X.mult(Xt);         // Gram matrix
        pMatrix alpha = G.add( G.I().mult(0.01) ).inv().mult(Y);   // Lagrange multipliers
        pMatrix w = Xt.mult(alpha); // weights
        */

        double[] r = new double[n+1];
        for(int i=0;i<r.length;i++)
            r[i] = w.get(i,0);
        return r;
    }

    /****************************************************************
     * IMAGE PROCESSING STUFF
     ****************************************************************/
    
    /**
     * convolve current matrix with kernel k.
     */
    public pMatrix convolve(pMatrix k){
        pMatrix v = new pMatrix(rows,cols);
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                double sum = 0;
                for(int ik=0;ik<k.rows;ik++){
                    int _i = (rows + i - k.rows/2 + ik) % rows;
                    for(int jk=0;jk<k.cols;jk++){
                        int _j = (cols + j - k.cols/2 + jk) % cols;
                        sum += arr[(_i << pitch) + _j] * k.arr[(ik << k.pitch) + jk];
                    }
                }
                v.arr[(i << v.pitch) + j] = sum;
            }
        }
        return v;
    }

    /****************************************************************
     * FACTORIZATION METHODS.
     ****************************************************************/

    /**
     * non-negative matrix factorization  非负矩阵分解
     * 
     * w.mult(h) = this
     * w is weights matrix.   w为权重矩阵。
     * h is feature matrix.   h为特征矩阵。
     *
     * k = size of factor matrices.  k 为因子矩阵的大小。
     * n is NUMBER OF ITERATIONS.    n为迭代次数。
     * e is ERROR (once we're below that, we can return)
     * 
     * returns w,h.
     */
    public pMatrix[] nmf(int k,int n,double e){
        pMatrix w = urandom(rows,k,1,n);
        pMatrix h = urandom(k,cols,1,n);
        for(int i=0;i<n;i++){
            // compute output.
            pMatrix wh = w.mult(h);
            double cost = distSq(wh);

            // if found solution
            if(cost < e)
                break;

            // update feature matrix.
            pMatrix wt = w.T();
            pMatrix hn = wt.mult(this);
            pMatrix hd = wt.mult(wh);
            h.smultEq(hn.smultEq(hd.sinvEq()));

            // update weights matrix
            pMatrix ht = h.T();
            pMatrix wn = mult(ht);
            pMatrix wd = w.mult(h).mult(ht);
            w.smultEq(wn.smultEq(wd.sinvEq()));
        }
        return new pMatrix[]{ w, h };
    }

    /****************************************************************
     * 3D, 4x4 MATRIX OPERATIONS (MOSTLY FOR GRAPHICS).
     ****************************************************************/

    /**
     * cross multiply.  交叉相乘
     */
    public pMatrix xmult3d(pMatrix v2){
        pMatrix v1 = this;
        pMatrix c = v1.dup();
        c.arr[0] = v1.arr[1]*v2.arr[2] - v1.arr[2]*v2.arr[1];
        c.arr[1] = v1.arr[2]*v2.arr[0] - v1.arr[0]*v2.arr[2];
        c.arr[2] = v1.arr[0]*v2.arr[1] - v1.arr[1]*v2.arr[0];
        return c;
    }

    /**
     * 3D, 4x4 matrix.
     * returns a rotation matrix around X axix.
     * Note: change sin term signs to change handedness.
     *
     * nx = x;
     * ny = x * 0 + y * cos + z * -sin;
     * nz = x * 0 + y * sin + z * cos;
     */
    public static pMatrix rotatex3d(double a){
        pMatrix tmp = I(4,4);
        double cos = (double)Math.cos(a);
        double sin = (double)Math.sin(a);
        tmp.arr[(1<<tmp.pitch)+1] = cos;
        tmp.arr[(1<<tmp.pitch)+2] = -sin;
        tmp.arr[(2<<tmp.pitch)+1] = sin;
        tmp.arr[(2<<tmp.pitch)+2] = cos;
        return tmp;
    }

    /**
     * 3D, 4x4 matrix.
     * returns a rotation matrix around Y axix.
     *
     * nx = x * cos + y * 0 + z * sin;
     * ny = y;
     * nz = x * -sin + y * 0 + z * cos;
     */
    public static pMatrix rotatey3d(double a){
        pMatrix tmp = I(4,4);
        double cos = (double)Math.cos(a);
        double sin = (double)Math.sin(a);
        tmp.arr[(0<<tmp.pitch)+0] = cos;
        tmp.arr[(0<<tmp.pitch)+2] = sin;
        tmp.arr[(2<<tmp.pitch)+0] = -sin;
        tmp.arr[(2<<tmp.pitch)+2] = cos;
        return tmp;
    }

    /**
     * 3D, 4x4 matrix.
     * returns a rotation matrix around Z axix.
     *
     * nx = x * cos + y * -sin + z * 0;
     * ny = x * sin + y * cos + z * 0;
     * nz = z;
     */
    public static pMatrix rotatez3d(double a){
        pMatrix tmp = I(4,4);
        double cos = (double)Math.cos(a);
        double sin = (double)Math.sin(a);
        tmp.arr[(0<<tmp.pitch)+0] = cos;
        tmp.arr[(0<<tmp.pitch)+1] = -sin;
        tmp.arr[(1<<tmp.pitch)+0] = sin;
        tmp.arr[(1<<tmp.pitch)+1] = cos;
        return tmp;
    }

    /**
     * 3D, 4x4 matrix.
     * returns a translation matrix 
     */
    public static pMatrix translate3d(double x,double y,double z){
        pMatrix tmp = I(4,4);
        tmp.arr[(0<<tmp.pitch)+3] = x;
        tmp.arr[(1<<tmp.pitch)+3] = y;
        tmp.arr[(2<<tmp.pitch)+3] = z;
        return tmp;
    }

    /**
     * 3D, 4x4 matrix.
     * returns a translation matrix 
     */
    public static pMatrix translate3d(double[] v){
        pMatrix tmp = I(4,4);
        tmp.arr[(0<<tmp.pitch)+3] = v[0];
        tmp.arr[(1<<tmp.pitch)+3] = v[1];
        tmp.arr[(2<<tmp.pitch)+3] = v[2];
        return tmp;
    }


    /**
     * 3D, 4x4 matrix.
     * scales the matrix.
     */
    public static pMatrix scale3d(double s){
        pMatrix tmp = I(4,4);
        tmp.arr[(0<<tmp.pitch)+0] = s;
        tmp.arr[(1<<tmp.pitch)+1] = s;
        tmp.arr[(2<<tmp.pitch)+2] = s;
        return tmp;
    }


    /**
     * 3D, 4x4 matrix.
     * scales the matrix.
     */
    public static pMatrix scale3d(double x,double y,double z){
        pMatrix tmp = I(4,4);
        tmp.arr[(0<<tmp.pitch)+0] = x;
        tmp.arr[(1<<tmp.pitch)+1] = y;
        tmp.arr[(2<<tmp.pitch)+2] = z;
        return tmp;
    }


    /**
     * 3D, 4x4 matrix.
     * scales a matrix in relation to a point
     *
     * @param s The scale.
     * @param v The center of point where to scale.
     */
    public static pMatrix scale3d(double s,double[] v){
        pMatrix tmp = I(4,4);
        tmp.arr[(0<<tmp.pitch)+0] = s;
        tmp.arr[(0<<tmp.pitch)+3] = (1 - s) * v[0];
        tmp.arr[(1<<tmp.pitch)+1] = s;
        tmp.arr[(1<<tmp.pitch)+3] = (1 - s) * v[1];
        tmp.arr[(2<<tmp.pitch)+2] = s;
        tmp.arr[(2<<tmp.pitch)+3] = (1 - s) * v[2];
        return tmp;
    }

    /**
     * 3D, 4x4 matrix.
     * scales a matrix in relation to a point
     *
     * @param s The scales (one for each coord).
     * @param v The center of point where to scale.
     */
    public static pMatrix scale3d(double[] s,double[] v){
        pMatrix tmp = I(4,4);
        tmp.arr[(0<<tmp.pitch)+0] = s[0];
        tmp.arr[(0<<tmp.pitch)+3] = (1 - s[0]) * v[0];
        tmp.arr[(1<<tmp.pitch)+1] = s[1];
        tmp.arr[(1<<tmp.pitch)+3] = (1 - s[1]) * v[1];
        tmp.arr[(2<<tmp.pitch)+2] = s[2];
        tmp.arr[(2<<tmp.pitch)+3] = (1 - s[2]) * v[2];
        return tmp;
    }


    /**
     * 3D, 4x4 matrix.    
     * reflect in the X axis.
     */
    public static pMatrix reflectx3d(){
        pMatrix tmp = I(4,4);
        tmp.arr[(0<<tmp.pitch)+0] = -tmp.arr[(0<<tmp.pitch)+0];
        return tmp;
    }

    /**
     * 3D, 4x4 matrix.    
     * reflect in the Y axis.
     */
    public static pMatrix reflecty3d(){
        pMatrix tmp = I(4,4);
        tmp.arr[(1<<tmp.pitch)+1] = -tmp.arr[(1<<tmp.pitch)+1];
        return tmp;
    }

    /**
     * 3D, 4x4 matrix.    
     * reflect in the Z axis.
     *
     * also known as conversion from right handled
     * to left handled coord system, or vice versa.
     */
    public static pMatrix reflectz3d(){
        pMatrix tmp = I(4,4);
        tmp.arr[(2<<tmp.pitch)+2] = -tmp.arr[(2<<tmp.pitch)+2];
        return tmp;
    }

    /**
     * 3D, 4x4 matrix.    
     * create perspective matrix; turns current matrix into
     * perspective matrix.
     * 
     * Note: changes current matrix!
     *
     * @param zrpr at what z is the eye?
     * @param zvp at what z is the view plane?
     */
    public pMatrix perspective3d(double zprp,double zvp){
        double t2,t3;
        double d1 = zprp - zvp;
        double a = -zvp/d1;
        double b = -1/d1;
        double c = zvp*(zprp/d1);
        double d = zvp/d1;
        for (int i=0;i<4;i++) {
            t2 = arr[(i<<pitch)+2];
            t3 = arr[(i<<pitch)+3];
            arr[(i<<pitch)+2] = t2*a + t3*b;
            arr[(i<<pitch)+3] = t2*c + t3*d;
        }
        return this;
    }

    /****************************************************************
     * QUATERNION STUFF [w,x,y,z] format.
     ****************************************************************/

    /**
     * axis angle form; x,y,z,theta
     */
    public static pMatrix toQuaternion(double ax,double ay,double az,double theta){
        double s = Math.sin(theta/2);
        return new pMatrix(new double[]{Math.cos(theta/2),s*ax,s*ay,s*az } );
    }

    /**
     * euler form; x-angle, y-angle, z-angle
     */
    public static pMatrix toQuaternion(double ax,double ay,double az){
        pMatrix qx,qy,qz,qt;
        qx = toQuaternion(1,0,0,ax);
        qy = toQuaternion(0,1,0,ay);
        qz = toQuaternion(0,0,1,az);
        return qx.quatMult(qy).quatMult(qz);
    }

    /**
     * Qr = Q1.Q2 = ( w1.w2 - v1.v2, w1.v2 + w2.v1 + v1 x v2 )
     */
    public pMatrix quatMult(pMatrix q){
        return quatMult(this,q);
    }

    /**
     * quaternion multiplication
     */
    public static pMatrix quatMult(pMatrix A,pMatrix B){
        int w=0,x=1,y=2,z=3;
        pMatrix C = new pMatrix(4,1);
        C.arr[x] = A.arr[w]*B.arr[x] + A.arr[x]*B.arr[w] + 
            A.arr[y]*B.arr[z] - A.arr[z]*B.arr[y];
        C.arr[y] = A.arr[w]*B.arr[y] - A.arr[x]*B.arr[z] + 
            A.arr[y]*B.arr[w] + A.arr[z]*B.arr[x];
        C.arr[z] = A.arr[w]*B.arr[z] + A.arr[x]*B.arr[y] - 
            A.arr[y]*B.arr[x] + A.arr[z]*B.arr[w];
        C.arr[w] = A.arr[w]*B.arr[w] - A.arr[x]*B.arr[x] - 
            A.arr[y]*B.arr[y] - A.arr[z]*B.arr[z];        
        return C;
    }

    /**
     * convert to 4x4 rotation matrix
     */
    public pMatrix quatToM(){
        double w = get(0), x=get(1), y=get(2), z=get(3);
        double[] m = {
            1 - 2*y*y - 2*z*z, 2*x*y - 2*w*z, 2*x*z + 2*w*y, 0,
            2*x*y + 2*w*z, 1 - 2*x*x-2*z*z, 2*y*z - 2*w*x, 0,
            2*x*z - 2*w*y,  2*y*z + 2*w*x,  1-2*x*x-2*y*y, 0,
            0, 0, 0, 1
        };
        pMatrix M = new pMatrix(4,4);
        System.arraycopy(m,0,M.arr,0,m.length);
        return M;
    }

    /**
     * convert to axis/angle form  转换为轴/角的形式
     * 
     * Ax, Ay, Az, Theta
     */
    public double[] quatToAxisAngle(){
        double w = get(0), x=get(1), y=get(2), z=get(3);        
        double s = x*x + y*y + z*z;
        double ax = x / s;
        double ay = y / s;
        double az = z / s;
        double theta = 2.0 * Math.acos(w);
        double[] aa = {ax,ay,az,theta};
        return aa;
    }

    /**
     * return the inverse
     */
    public pMatrix quatConjugate(){
        pMatrix tmp = dup();
        for(int i=1;i<=3;i++)
            tmp.arr[i] = -tmp.arr[i];
        return tmp;
    }

    /****************************************************************
     * 3D MATRIX INVERSE STUFF.
     ****************************************************************/

    /**
     * 3D, 4x4 matrix.    
     * returns the inverse of the current matrix
     */
    public pMatrix inv3d(){
        pMatrix n = I(4,4);
        if(m4_inverse(n.arr,arr)){
            return n;
        }
        return I(4,4);
    }

    /**
     * various matrix operations 
     * (inverse; probably faster than LU decomposition 
     * method above.). taken from:
     * http://skal.planet-d.net/demo/matrixfaq.htm
     */

    /**
     * compute determinant of a 3x3 matrix
     */
    private static double m3_det(double[] mat){
        double det;
        det = mat[0] * ( mat[4]*mat[8] - mat[7]*mat[5] )
                - mat[1] * ( mat[3]*mat[8] - mat[6]*mat[5] )
                + mat[2] * ( mat[3]*mat[7] - mat[6]*mat[4] );
        return det;
    }

    private static boolean m3_inverse(double[] mr,double[] ma){
        double det = m3_det( ma );
        if ( Math.abs( det ) < (double)0.0005 ) {
            //m3_identity( mr );
            return false;
        }
        mr[0] =    ma[4]*ma[8] - ma[5]*ma[7]   / det;
        mr[1] = -( ma[1]*ma[8] - ma[7]*ma[2] ) / det;
        mr[2] =    ma[1]*ma[5] - ma[4]*ma[2]   / det;
        mr[3] = -( ma[3]*ma[8] - ma[5]*ma[6] ) / det;
        mr[4] =    ma[0]*ma[8] - ma[6]*ma[2]   / det;
        mr[5] = -( ma[0]*ma[5] - ma[3]*ma[2] ) / det;
        mr[6] =    ma[3]*ma[7] - ma[6]*ma[4]   / det;
        mr[7] = -( ma[0]*ma[7] - ma[6]*ma[1] ) / det;
        mr[8] =    ma[0]*ma[4] - ma[1]*ma[3]   / det;
        return true;
    }

    /**
     * mr is 4x4, mb is 3x4
     */
    private static void m4_submat(double[] mr,double[] mb, int i, int j ) {
        int di, dj, si, sj;
        // loop through 3x3 submatrix
        for( di = 0; di < 3; di ++ ) {
            for( dj = 0; dj < 3; dj ++ ) {
                // map 3x3 element (destination) to 4x4 element (source)
                si = di + ( ( di >= i ) ? 1 : 0 );
                sj = dj + ( ( dj >= j ) ? 1 : 0 );
                // copy element
                mb[di * 3 + dj] = mr[si * 4 + sj];
            }
        }
    }

    private static double m4_det(double[] mr) {
        double det, result = 0, i = 1;
        double[] msub3 = new double[3*3];
        int     n;
        for ( n = 0; n < 4; n++, i *= -1 ){
            m4_submat( mr, msub3, 0, n );
            det     = m3_det( msub3 );
            result += mr[n] * det * i;
        }
        return result;
    }

    private static boolean m4_inverse(double[] mr, double[] ma ){
        double mdet = m4_det( ma );
        double[] mtemp = new double[3*3];
        int     i, j, sign;
        if (Math.abs( mdet ) < (double)0.0005 ){
            // m4_identity( mr );
            return false;
        }
        for ( i = 0; i < 4; i++ )
            for ( j = 0; j < 4; j++ ){
                sign = 1 - ( (i +j) % 2 ) * 2;
                m4_submat( ma, mtemp, i, j );
                mr[i+j*4] = ( m3_det( mtemp ) * sign ) / mdet;
            }
        return true;
    }

    
    /****************************************************************
     * JAVA IMAGE RELATED BLOAT
     ****************************************************************/

    /**
     * make an image producer; so someone can easily turn this into a
     * gray scale image, ie:
     * g.drawImage(createImage(m.getImageProducer()),0,0,null);
     */

    public java.awt.image.ImageProducer getImageProducer(){
        return getImageProducer(new pMatrix[]{this,this,this});
    }

    /**
     * given any 3 matrices (of same size) combine them via rgb into single image.
     * Note that these are static images (not animating).
     */
    public static java.awt.image.ImageProducer getImageProducer(final pMatrix[] rgb){
        return new java.awt.image.ImageProducer(){
            java.util.Vector v = new java.util.Vector();

            public void addConsumer(java.awt.image.ImageConsumer ic){
               if(!isConsumer(ic))
                    v.addElement(ic);
              requestTopDownLeftRightResend(ic);
            }
           public boolean isConsumer(java.awt.image.ImageConsumer ic){
                return v.contains(ic);
            }
            public void removeConsumer(java.awt.image.ImageConsumer ic){
                v.removeElement(ic);
            }
           public void startProduction(java.awt.image.ImageConsumer ic){
                addConsumer(ic);
            }
            public void requestTopDownLeftRightResend(java.awt.image.ImageConsumer ic){
               int[] pixels = new int[rgb[0].arr.length];
               java.awt.image.ColorModel cm = java.awt.image.ColorModel.getRGBdefault();
                ic.setDimensions(rgb[0].cols,rgb[0].rows);
                ic.setColorModel(cm);
                ic.setHints(
                    java.awt.image.ImageConsumer.TOPDOWNLEFTRIGHT |
                    java.awt.image.ImageConsumer.COMPLETESCANLINES |
                    java.awt.image.ImageConsumer.SINGLEPASS
               );
                for(int i=0;i<rgb[0].rows;i++){
                    for(int j=0;j<rgb[0].cols;j++){
                        pixels[(i<<rgb[0].pitch)+j] = 0xFF << 24 | 
                            ((int)(rgb[0].arr[(i<<rgb[0].pitch)+j]) % 0x100) << 16 |                            ((int)(rgb[1].arr[(i<<rgb[1].pitch)+j]) % 0x100) << 8 | 
                            ((int)(rgb[2].arr[(i<<rgb[2].pitch)+j]) % 0x100);
                    }
                }
               ic.setPixels(0,0,rgb[0].cols,rgb[0].rows,cm,pixels,0,1<<rgb[0].pitch);
                ic.imageComplete(java.awt.image.ImageConsumer.STATICIMAGEDONE);
            }
        };
    }
    
    /**
     * get pixels from an image; RBG
     */
    /*
    public static pMatrix[] getRGB(java.awt.Image img){
        

    MediaTracker tracker;
    
    // double buffer stuff.
    int W, H, pix[];
    Image im;
    MemoryImageSource mis;
        // load the image
        image = getImage(getDocumentBase(),param);
        if(image == null)
            return;
            
        tracker = new MediaTracker(this);
        tracker.addImage(image, 0);
        
        // grab pixels.
        try {
            tracker.waitForID(0);
            m_width = image.getWidth(this);
            m_height = image.getHeight(this);
            pixels = new int[m_width * m_height];
            PixelGrabber pg = new PixelGrabber(image,0,0,m_width,m_height,pixels,0,m_width);
            pg.grabPixels();            
        } catch (InterruptedException e) {
            return;
        }
        
    }
    */

    /****************************************************************
     * TO STRING.
     ****************************************************************/

    /**
     * output matrix in csv format
     */
    public String toString(){
        StringBuffer sb = new StringBuffer();
        sb.append("pMatrix["+rows+"x"+cols+"][pitch="+pitch+"][len="+arr.length+"]:\n");
        /*
        NumberFormat nf = NumberFormat.getInstance();
        nf.setGroupingUsed(false);
        nf.setMaximumFractionDigits(0x2);
        nf.setMaximumIntegerDigits(0xFF);
        nf.setMinimumFractionDigits(0);
        nf.setMinimumIntegerDigits(1);
        */

        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                //sb.append(nf.format(get(i,j)) + (j+1 < cols ? "," : "\n") );
                sb.append(get(i,j) + (j+1 < cols ? "," : "\n") );
            }
        }
        return sb.toString();
    }



    /****************************************************************
     * DIRTY TEST CODE
     ****************************************************************/

    // quickly test stuff.
    public static void main(String[] args) throws Exception {
        //pMatrix a = pMatrix.nrandom(5,5,50,25);
    	double test1[] = new double[]{22,28};
    	double test2[] = new double[]{22,28};
        pMatrix a = new pMatrix(new double[][]{
           test1,
           test2,
           test1,
           test2
        });

        pMatrix[] b = a.nmf(3,100,0);
        System.out.println("a: "+a);
        System.out.println("b[0]: "+b[0]);
        System.out.println("b[1]: "+b[1]);
        pMatrix p = b[0].mult(b[1]);
        System.out.println("b[0].mult(b[1]): "+p);
        double dist = a.dist(p);
        System.out.println("dist: "+dist);
        
    }
}


