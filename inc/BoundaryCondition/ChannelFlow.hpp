
// *  Channel flow
// *  Velocities: All Dirichlet (all 0 except u = u0 only at West) 
// *  except outlet (East).
// *  Pressure:   All Neumann.

    Dir[0][0] = 1.;/*   */Num[0][0] = 0.;
    Dir[0][1] = 0.;/* 0 */Num[0][1] = 0.;
    Dir[0][2] = 0.;/*   */Num[0][2] = 0.;

    Dir[1][0] = 0.;/*   */Num[1][0] = 1.;
    Dir[1][1] = 0.;/* 1 */Num[1][1] = 1.;
    Dir[1][2] = 0.;/*   */Num[1][2] = 1.;


    Dir[2][0] = 0.;/*   */Num[2][0] = 0.;
    Dir[2][1] = 0.;/* 2 */Num[2][1] = 0.;
    Dir[2][2] = 0.;/*   */Num[2][2] = 0.;

    Dir[3][0] = 0.;/*   */Num[3][0] = 0.;
    Dir[3][1] = 0.;/* 3 */Num[3][1] = 0.;
    Dir[3][2] = 0.;/*   */Num[3][2] = 0.;


    Dir[4][0] = 0.;/*   */Num[4][0] = 0.;
    Dir[4][1] = 0.;/* 4 */Num[4][1] = 0.;
    Dir[4][2] = 0.;/*   */Num[4][2] = 0.;


    Dir[5][0] = 0.;/*   */Num[5][0] = 0.;
    Dir[5][1] = 0.;/* 5 */Num[5][1] = 0.;
    Dir[5][2] = 0.;/*   */Num[5][2] = 0.;
