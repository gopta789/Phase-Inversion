#include <bits/stdc++.h>
//#include"Elasticity.h"
using namespace std;

/* ================= PARAMETERS ================= */
const int nx = 128;
const int ny = 128;

const double dx = 1.0;
const double dt = 0.005;

const double Lmob = 1.0;
const double kappa_g = 1.0;
const double P = 1.0;  

const double eT      = 0.01;    // eigenstrain          [Table 1, both papers]
const double DX     = 1.0;     // grid spacing x       [Table 1, ref2]
const double DY      = 1.0;     // grid spacing y       [Table 1, ref2]
const double TOL     = 1e-8;    // convergence tol      [ref2 Eq.43]
const int MAX_ITER   = 2;  


/* ================= FIELDS ================= */
double g[nx][ny];
double g_new[nx][ny];

double UX[nx][ny], UY[nx][ny];
//double g[nx][ny];

/* stiffness tensors from paper */
double Cm[2][2][2][2];
double Cp[2][2][2][2];
double Ceff[2][2][2][2];
double DeltaC[2][2][2][2];


/* transformation strain */
double epsT=0.01;

void init_composition(){
    const double R0    = 20.0;   // precipitate radius
    const double width = 3.0;    // diffuse interface half-width (grid units)
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double r = sqrt((i - nx/2.0)*(i - nx/2.0) + (j - ny/2.0)*(j - ny/2.0));
            // tanh profile: g=1 inside, g=0 outside, smooth transition at r=R0
            g[i][j] = 0.5 * (1.0 - tanh((r - R0) / width));
        }
    }
}
double b_func(double cv)
{
    return cv*cv*cv * (10.0 - 15.0*cv + 6.0*cv*cv);
}

void init_Ceff()
{
    // --- zero out all ---
    for(int i=0;i<2;i++)
    for(int j=0;j<2;j++)
    for(int k=0;k<2;k++)
    for(int l=0;l<2;l++){
        Cm[i][j][k][l]   = 0.0;
        Cp[i][j][k][l]   = 0.0;
        Ceff[i][j][k][l] = 0.0;
    }

    // --- C^m components ---
    // Replace with actual values from your material system
    Cm[0][0][0][0] = 1.0;    // C^m_1111
    Cm[0][0][1][1] = 0.5;    // C^m_1122
    Cm[1][1][0][0] = 0.5;    // C^m_2211
    Cm[1][1][1][1] = 1.0;    // C^m_2222
    Cm[0][1][0][1] = 0.3;    // C^m_1212
    Cm[0][1][1][0] = 0.3;    // C^m_1221
    Cm[1][0][0][1] = 0.3;    // C^m_2112
    Cm[1][0][1][0] = 0.3;    // C^m_2121

    // --- C^p components for delta = Gp/Gm = 0.5 (soft precipitate) ---
    // [ref1 & ref2, delta = 0.5 case]
    Cp[0][0][0][0] = 0.5;    // C^p_1111
    Cp[0][0][1][1] = 0.25;   // C^p_1122
    Cp[1][1][0][0] = 0.25;   // C^p_2211
    Cp[1][1][1][1] = 0.5;    // C^p_2222
    Cp[0][1][0][1] = 0.15;   // C^p_1212
    Cp[0][1][1][0] = 0.15;   // C^p_1221
    Cp[1][0][0][1] = 0.15;   // C^p_2112
    Cp[1][0][1][0] = 0.15;   // C^p_2121

    // --- C^eff = 0.5*(C^m + C^p) ---
    // [Table 1, both papers]
    for(int i=0;i<2;i++)
    for(int j=0;j<2;j++)
    for(int k=0;k<2;k++)
    for(int l=0;l<2;l++){
        Ceff[i][j][k][l]    = 0.5*(Cm[i][j][k][l] + Cp[i][j][k][l]);
        DeltaC[i][j][k][l]  = Cp[i][j][k][l] - Cm[i][j][k][l];  // BUG FIX: was never set
    }
}

void init_displacement()
{
    for(int i=0;i<nx;i++)
    for(int j=0;j<ny;j++){
        UX[i][j] = 0.0;
        UY[i][j] = 0.0;
    }
}

double getU(int l, int ii, int jj)
{
    if(ii <= 0 || ii >= nx-1 || jj <= 0 || jj >= ny-1)
        return 0.0;
    return (l == 0) ? UX[ii][jj] : UY[ii][jj];
}

double diag_coeff(int i_eq, int l)
{
    double d = 0.0;
    d += Ceff[i_eq][0][0][l] * (-2.0 / (DX*DX));
    d += Ceff[i_eq][1][1][l] * (-2.0 / (DY*DY));
    return d;
}

double offdiag_LHS(int i_eq, int p, int q)
{
    double val = 0.0;

    for(int j=0; j<2; j++)
    for(int k=0; k<2; k++)
    for(int l=0; l<2; l++)
    {
        double coeff = Ceff[i_eq][j][k][l];
        if(coeff == 0.0) continue;

        if(j==0 && k==0)
        {
            // d^2 Ul/dr1^2 — neighbor terms only, exclude -2*Ul(p,q)
            val += coeff * ( getU(l,p+1,q) + getU(l,p-1,q) )
                         / (DX*DX);
        }
        else if(j==1 && k==1)
        {
            // d^2 Ul/dr2^2 — neighbor terms only, exclude -2*Ul(p,q)
            val += coeff * ( getU(l,p,q+1) + getU(l,p,q-1) )
                         / (DY*DY);
        }
        else if((j==0 && k==1) || (j==1 && k==0))
        {
            // d^2 Ul/dr1 dr2 — all 4 corners, no center point
            val += coeff * ( getU(l,p+1,q+1) - getU(l,p+1,q-1)
                            -getU(l,p-1,q+1) + getU(l,p-1,q-1) )
                         / (4.0*DX*DY);
        }
    }

    return val;
}

double compute_RHS(int i_eq, int p, int q)
{
    // Safe getter for composition at boundary
    auto getC = [&](int ii, int jj) -> double {
        int ic = max(0, min(nx-1, ii));
        int jc = max(0, min(ny-1, jj));
        return g[ic][jc];
    };

    // db/dr_1 : central difference
    double dbdr1 = (b_func(getC(p+1,q)) - b_func(getC(p-1,q)))
                   / (2.0*DX);

    // db/dr_2 : central difference
    double dbdr2 = (b_func(getC(p,q+1)) - b_func(getC(p,q-1)))
                   / (2.0*DY);

    // Contraction: Ceff2[i][j][k][k] summed over k
    // j=0 (r_1 direction): Ceff2[i][0][0][0] + Ceff2[i][0][1][1]
    // j=1 (r_2 direction): Ceff2[i][1][0][0] + Ceff2[i][1][1][1]
    double coeff_j0 = Ceff[i_eq][0][0][0] + Ceff[i_eq][0][1][1];
    double coeff_j1 = Ceff[i_eq][1][0][0] + Ceff[i_eq][1][1][1];

    return eT * (coeff_j0 * dbdr1 + coeff_j1 * dbdr2);
}

double gauss_seidel_sweep()
{
    double error = 0.0;

    for(int p=1; p<nx-1; p++)
    for(int q=1; q<ny-1; q++)
    {
        // --- Diagonal coefficients of 2x2 system ---
        double D00 = diag_coeff(0, 0);   // coeff of UX in eq i=0
        double D01 = diag_coeff(0, 1);   // coeff of UY in eq i=0
        double D10 = diag_coeff(1, 0);   // coeff of UX in eq i=1
        double D11 = diag_coeff(1, 1);   // coeff of UY in eq i=1

        // --- Off-diagonal neighbor contributions ---
        double off0 = offdiag_LHS(0, p, q);
        double off1 = offdiag_LHS(1, p, q);

        // --- RHS ---
        double R0 = compute_RHS(0, p, q);
        double R1 = compute_RHS(1, p, q);

        // --- Modified RHS: move neighbor terms to RHS ---
        double b0 = R0 - off0;
        double b1 = R1 - off1;

        // --- Solve 2x2 by Cramer's rule ---
        double det = D00*D11 - D01*D10;

        double UX_old = UX[p][q];
        double UY_old = UY[p][q];

        if(abs(det) > 1e-15)
        {
            UX[p][q] = (b0*D11 - b1*D01) / det;
            UY[p][q] = (D00*b1 - D10*b0) / det;
        }

        // --- Accumulate L2 error [ref2 Eq.43] ---
        double dUX = UX[p][q] - UX_old;
        double dUY = UY[p][q] - UY_old;
        error += dUX*dUX + dUY*dUY;
    }

    return sqrt(error);
}
void solve_elasticity()
{
   
    double error = 1.0;
    int iter = 0;

    while(error > TOL && iter < MAX_ITER)
    {
        error = gauss_seidel_sweep();
        iter++;
    }
  
}

double beta(double x)
{
    return x*x*x*(10-15*x+6*x*x);
}

double phi(double x)
{
    return beta(x);   // same interpolation used in paper
}
double beta_prime(double x){
    if (x <= 0.0 || x >= 1.0) return 0.0;
    return 30.0*x*x*(1.0 - x)*(1.0 - x);
}
double phi_prime(double x){
    if (x <= 0.0 || x >= 1.0) return 0.0;
    return 30.0*x*x*(1.0 - x)*(1.0 - x);
}

/* interpolation */

/* ================= FUNCTIONS ================= */

// W(g)
double W(double g) {
    if (g <= 0.0) return 0.0;
    if (g >= 1.0) return 1.0;
    return g*g*g*(10.0 - 15.0*g + 6.0*g*g);
}

// W'(g)
double dWdg(double g) {
    if (g <= 0.0 || g >= 1.0) return 0.0;
    return 30.0*g*g*(1.0 - g)*(1.0 - g);
}


double laplacian(int i, int j) {
    //cout<<i<<" "<<j<<endl;

    int ip = (i + 1);
    int im = (i - 1);
    int jp = (j + 1);
    int jm = (j - 1);

    if(i==0){
        im = (i);
    }
    if(i==nx-1){
        ip = (i);
    }
    if(j==0){
        jm = (j);
    }
    if(j==ny-1){
        jp = (j);
    }  

    return ( (
        g[ip][j] + g[im][j] +
        g[i][jp] + g[i][jm] -
        4.0 * g[i][j]
    ) / (dx * dx) );
}
/* ================= Elastic Part start ================= */

/* fields */

/* compute stress using Eq 13 + Eq 16 */
void compute_stress(
double exx_el,double eyy_el,double exy_el,
double gval,
double &sxx,double &syy,double &sxy)
{
    double eps[2][2];

    eps[0][0]=exx_el;
    eps[1][1]=eyy_el;
    eps[0][1]=exy_el;
    eps[1][0]=exy_el;

    double sigma[2][2]={0};

    for(int i=0;i<2;i++)
    for(int j=0;j<2;j++)
    for(int k=0;k<2;k++)
    for(int l=0;l<2;l++)
    {
        double Cijkl =
        Ceff[i][j][k][l]
        +phi(gval)*DeltaC[i][j][k][l];

        sigma[i][j]+=Cijkl*eps[k][l];
    }

    sxx=sigma[0][0];
    syy=sigma[1][1];
    sxy=sigma[0][1];
}

/* compute elastic energy Eq 10 */
double compute_elastic_energy(int i,int j)
{
    // BUG FIX: central-difference stencil needs one ghost layer on each side
    // — skip boundary cells to avoid out-of-bounds UX/UY access (was giving NaN)
    if(i==0 || i==nx-1 || j==0 || j==ny-1) return 0.0;

    double Fel=0;

        /* strain Eq 12 */

    double exx=(UX[i+1][j]-UX[i-1][j])/(2*dx);
    double eyy=(UY[i][j+1]-UY[i][j-1])/(2*dx);

    double exy=0.5*((UX[i][j+1]-UX[i][j-1])/(2*dx)
                       +(UY[i+1][j]-UY[i-1][j])/(2*dx));

        /* eigenstrain Eq 15 */

    // cout<<"exx: "<<exx<<endl;
    // cout<<"eyy: "<<eyy<<endl;
    // cout<<"exy: "<<exy<<endl;

    double e0=beta(g[i][j])*epsT;

    double exx_el=exx-e0;
    double eyy_el=eyy-e0;
    double exy_el=exy;

        /* stress Eq 13 + Eq 16 */

    double sxx,syy,sxy;
    compute_stress(exx_el,eyy_el,exy_el,g[i][j],sxx,syy,sxy);

        /* energy density Eq 10 */
    /* ---------- stress trace ---------- */

    double sigma_trace = sxx + syy;

    /* ---------- tensor contraction term ---------- */

    double strain_contract = 0.0;

    double eps_el[2][2];

    eps_el[0][0]=exx_el;
    eps_el[1][1]=eyy_el;
    eps_el[0][1]=exy_el;
    eps_el[1][0]=exy_el;

    for(int a=0;a<2;a++)
    for(int b=0;b<2;b++)
    for(int c=0;c<2;c++)
    for(int d=0;d<2;d++)
    {
        strain_contract +=
        DeltaC[a][b][c][d] *
        eps_el[a][b] *
        eps_el[c][d];
    }

    /* ---------- final elastic term ---------- */

    double term =
        0.5 * (phi_prime(g[i][j]) * strain_contract
        - beta_prime(g[i][j]) * epsT * sigma_trace);

    cout<<"Term: "<<term<<endl;
    return term;
}

/* ================= Elastic Part End ================= */

/* ================= MAIN ================= */
int main() {

    init_Ceff();
    init_composition();
    init_displacement();
    ofstream fout_init("C:/Users/Aashu/plotting/eq2_initial_elastic.dat");
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            fout_init << i << " " << j << " " << g[i][j] << "\n";
        }
        fout_init << "\n";
    }
    fout_init.close();

    const double Achem = -1.0;

    // Time loop
    double el=0.0;
    for (int step = 0; step < 5; step++) {
        solve_elasticity(); //Elastic Part
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) { 

                double gij = g[i][j];

                // Chemical driving force and double well
                double chem =
                    Achem * dWdg(gij)
                    + 2.0 * P * gij * (1.0 - gij) * (1.0 - 2.0 * gij);
                //Elastic Part
                el=compute_elastic_energy(i,j); //elastic_AC(i,j);

                // Gradient term
                double grad = -2.0 * kappa_g * laplacian(i, j);

                // Allen–Cahn update
                g_new[i][j] =
                    gij - dt * Lmob * (chem + grad + el);
            }
        }

        // Update field
        //memcpy(g, g_new, sizeof(g));
        for(int i=0;i<nx;++i){
            for(int j=0;j<ny;++j){
                g[i][j]=g_new[i][j];
            }
        }

        // Output occasionally
        cout<<"Step"<<step<<endl;
        // if (step % 500 == 0) {
        //     cout << "Step " << step << endl;
        // }
        // if(step % 1000 == 0){
        //     char fname[100];
        //     sprintf(fname,"allen_results/g_%05d.dat",step);

        //     ofstream fout(fname);
        //     for(int i=0;i<nx;i++){
        //         for(int j=0;j<ny;j++)
        //         fout<<i<<" "<<j<<" "<<g[i][j]<<"\n";
        //         fout<<"\n";
        //     }
        //     fout.close();
        // }

    }
    ofstream fout_fnl("C:/Users/Aashu/plotting/eq2_final_elastic.dat");
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            fout_fnl << i << " " << j << " " << g[i][j] << "\n";
        }
        fout_fnl << "\n";
    }
    fout_fnl.close();

    // int ic = nx/2;
    // int jc = ny/2;

    // /* detect particle radius (g=0.5 crossing) */
    // double R=1.0;
    // for(int i=ic;i<nx-1;i++)
    // {
    //     if(g[i][jc]>=0.5 && g[i+1][jc]<0.5)
    //     {
    //         R = fabs(i-ic);
    //         break;
    //     }
    // }

    // ofstream fout_graph1("fig3_profile.dat");

    // for(int i=0;i<nx;i++)
    // {
    //     double x = i-ic;
    //     double xs = x/R;   // scaled coordinate
    //     fout<<xs<<" "<<g[i][jc]<<"\n";
    // }
    // fout_graph1.close();
    // cout<<"Fig3 profile generated"<<endl;

    return 0;
}
