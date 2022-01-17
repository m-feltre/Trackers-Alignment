

 #include "Math/Minimizer.h"
 #include "Math/Factory.h"
 #include "Math/Functor.h"
 #include "TRandom3.h"
 #include "TError.h"
 #include <iostream>
 #include <fstream>
 #include "TH1F.h"
 #include "TH2D.h"
 #include "TF1.h"
 #include "TCanvas.h"

TRandom3 *r=new TRandom3(42);
//INSTRUCTIONS: put in the same directory of "gaussian_output_hit_y_dy.txt" and "gaussian_output_lines_y_dy.txt"
//Create a folder named as "namefile" at line 48
//compile and call main() for minimization and print_all_it0() for plots


/* double displacement1=r->Uniform(-0.5,0.5);
double displacement2=r->Uniform(-0.5,0.5);
double displacement3=r->Uniform(-0.5,0.5);

double displ1=r->Uniform(-0.5,0.5);
double displ2=r->Uniform(-0.5,0.5);
double displ3=r->Uniform(-0.5,0.5); */
//Define the displacement for initializing the values for minimization. ->If all =0 they are at their true value
//displacement for x2, x3, x4 
double displacement1=0.2;
double displacement2=0.4;
double displacement3=0;
//displacement for y2, y3, y4
double displ1=0;
double displ2=0;
double displ3=0;
//Initializing variables m, q and a vector with 4 hits in y. 
   Double_t q=-3.82249;
   Double_t m=0.024043;
   Double_t y[4]={-3.82249, -2.8752, 2.97687, 4.03476 };

/////////////////////////////////////////////////////////Number of tracks per fit
   int n=8; //number of data processed. n=4,8 or 16. For 16 it is necessary to go on line 536 and change its value to be less than 2000, otherwise use 4000 max
    //Indeed the number of tracks is 32000->2000*16 or 4000*8 at MAXIMUM

   int contatore_tracce=0;

   std::string namefile="minimizer/gaussian_all_y/";
   //std::string namefile="minimizer/gaussian_x_lines/";


///////////////////////////////////////////True values
    Float_t x[4]={0.,39.4,282.8,326.8};
    Float_t dy[4]={0.,0.1,-0.3,0.9};
///////////////////////////////////////////Values after iteration 0
    Float_t result_iter_0[]={0.,0.,0.,0.,0.,0.};
    std::vector<double> x2_it0, x3_it0, x4_it0, y2_it0, y3_it0, y4_it0;



    //double w_tot=39.4+282.8+326.8;
    //double w_tot=39.4*39.4+282.8*282.6+326.8*326.8;
    //double w2=39.4, w3=282.8, w4=326.8;
    //double w2=39.4*39.4, w3=282.8*282.6, w4=326.8*326.8;

    double w_tot=1.;
    double w2=1., w3=1., w4=1.;


    ///////A struct dato contains 4 hit (1 per tracker) and the parameters of the line

    struct dato {
    std::vector<double> y={0.,0.,0.,0.};
    double m=0;
    double q=0;
    };

    std::vector<dato> dataset;

    //Mean of a vector
 double mean_v(std::vector<double> v) {
     double sum=0;
     for(int i=0;i<v.size();i++) {
         sum+=v.at(i);
     }
     sum=sum/v.size();
     return sum;
 }
 
 //////////////////////////////////////Function to be minimized
 double chi_square(const double *xx )
 {
    


   Double_t chisq = 0;
   Double_t delta=0;
   Double_t sum1=0, sum2=0, sum3=0, sum4=0;
   const Double_t x2 = xx[0];
   const Double_t x3 = xx[1];
   const Double_t x4 = xx[2];

   const Double_t dy2 = xx[3];
   const Double_t dy3 = xx[4];
   const Double_t dy4 = xx[5];

    /////// loop "for" on a dataset: its size depends on n (line 42)-> it is recreated for every "experiment" in the main function
   for(int i=0; i<dataset.size();i++) {
       dato a=dataset.at(i);
       y[0]=a.y[0];
       y[1]=a.y[1];
       y[2]=a.y[2];
       y[3]=a.y[3];
       /*
       m=a.m;
       q=a.q;
       */
     
      const Double_t m=xx[5+2*i+1];
      const Double_t q=xx[5+2*i+2];

        ////in this case w e w tot are equal to 1.->they are not relevant

         sum1=(y[0]-q)*(y[0]-q)/(1*1);
         sum2=w2/w_tot * (y[1]-q-m*x2+dy2)*(y[1]-q-m*x2+dy2);
         sum3=w3/w_tot*(y[2]-q-m*x3+dy3)*(y[2]-q-m*x3+dy3);
         sum4=w4/w_tot*(y[3]-q-m*x4+dy4)*(y[3]-q-m*x4+dy4); 

         //std::cerr<<sum1<<" "<<sum2<<" "<<sum3<<" "<<sum4<<std::endl;
   
        delta=delta+sum1+sum2+sum3+sum4;

   }

   
    //I found that for a given power the search for the minimum is more quick
    return pow(delta,2);
 }
 


 ///////Minimizer algorithm:
 /*
 It repeats two times:
    iteration=0
->First search of the minimum values

    iteration=1
->Second search in which the values of initialization are the mean for x2, x3, x4 and y2, y3, y4 found in iteration=0
->It does not seem to work 

 */
 std::vector<double> NumericalMinimization(int iteration=0)
 {
     std::vector<double> res_iteration={0.,0.,0.,0.,0.,0.};
if(iteration==0) {
    contatore_tracce++;
    
    //ROOT::Math::Minimizer* minimum =
    //   ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
 
    ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "");

    //I tried both with Minuit2 and GLSMultimin->The second seems better but I don't know why

    // set tolerance , etc...
    minimum->SetMaxFunctionCalls(100000000); // for Minuit/Minuit2
    minimum->SetMaxIterations(10000);  // for GSL
    minimum->SetTolerance(0.00000000001);
    minimum->SetPrintLevel(1);
    minimum->SetStrategy(2);
  
 
    // create function wrapper for minimizer
    // a IMultiGenFunction type
    ROOT::Math::Functor f(&chi_square,6+2*n);
    double step[5] = {0.01, 0.01, 0.01,0.0001, 0.00001};
    // starting point
 
    //double variable[3] = { 39.4, 282.8, 326.8};

    double variable[3] = { x[1]+displacement1,x[2]+displacement2,x[3]+displacement3};
    double lines[16] = { 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,};

 
    minimum->SetFunction(f);
 
    // Set the free variables to be minimized !
    //x positions
    minimum->SetVariable(0,"x2",variable[0], step[0]);
    minimum->SetVariable(1,"x3",variable[1], step[1]);
    minimum->SetVariable(2,"x4",variable[2], step[2]);

    minimum->SetVariableLimits(0,x[1]-2.,x[1]+2.);
    minimum->SetVariableLimits(1,x[2]-2.,x[2]+2.);
    minimum->SetVariableLimits(2,x[3]-2.,x[3]+2.);

    //dy of trackers

    

    minimum->SetVariable(3,"dy2",displ1+0.1, step[3]);
    minimum->SetVariable(4,"dy3",displ2-0.3, step[3]);
    minimum->SetVariable(5,"dy4",displ3+0.9, step[3]);

    minimum->SetVariableLimits(3,-0.9,+1.1);
    minimum->SetVariableLimits(4,-1.3,0.7);
    minimum->SetVariableLimits(5,-0.1,+1.9);

    
    //lines

    minimum->SetVariable(6,"m1",lines[1], step[4]);
    minimum->SetVariable(7,"q1",lines[1], step[0]);

    minimum->SetVariableLimits(6,-0.03,0.03);
    minimum->SetVariableLimits(7,-5.,5.);
    
    if(n>1) {
        minimum->SetVariable(8,"m2",lines[1], step[4]);
        minimum->SetVariable(9,"q2",lines[1], step[0]);

        minimum->SetVariableLimits(8,-0.03,0.03);
        minimum->SetVariableLimits(9,-5.,5.);
    }
    if(n>2) { //For 4 lines. It takes these two lines and all the previous one
        minimum->SetVariable(10,"m3",lines[1], step[4]);
        minimum->SetVariable(11,"q3",lines[1], step[0]);
        minimum->SetVariable(12,"m4",lines[1], step[4]);
        minimum->SetVariable(13,"q4",lines[1], step[0]);

        minimum->SetVariableLimits(10,-0.03,0.03);
        minimum->SetVariableLimits(11,-5.,5.);
        minimum->SetVariableLimits(12,-0.03,0.03);
        minimum->SetVariableLimits(13,-5.,5.);
    }
    if(n>4) {
        minimum->SetVariable(14,"m5",lines[1], step[4]);
        minimum->SetVariable(15,"q5",lines[1], step[0]);
        minimum->SetVariable(16,"m6",lines[1], step[4]);
        minimum->SetVariable(17,"q6",lines[1], step[0]);
        minimum->SetVariable(18,"m7",lines[1], step[4]);
        minimum->SetVariable(19,"q7",lines[1], step[0]);
        minimum->SetVariable(20,"m8",lines[1], step[4]);
        minimum->SetVariable(21,"q8",lines[1], step[0]);

        minimum->SetVariableLimits(14,-0.03,0.03);
        minimum->SetVariableLimits(15,-5.,5.);
        minimum->SetVariableLimits(16,-0.03,0.03);
        minimum->SetVariableLimits(17,-5.,5.);
        minimum->SetVariableLimits(18,-0.03,0.03);
        minimum->SetVariableLimits(19,-5.,5.);
        minimum->SetVariableLimits(20,-0.03,0.03);
        minimum->SetVariableLimits(21,-5.,5.);
    }
    if(n>8) {
        minimum->SetVariable(22,"m9",lines[1], step[4]);
        minimum->SetVariable(23,"q9",lines[1], step[0]);
        minimum->SetVariable(24,"m10",lines[1], step[4]);
        minimum->SetVariable(25,"q10",lines[1], step[0]);
        minimum->SetVariable(26,"m11",lines[1], step[4]);
        minimum->SetVariable(27,"q11",lines[1], step[0]);
        minimum->SetVariable(28,"m12",lines[1], step[4]);
        minimum->SetVariable(29,"q12",lines[1], step[0]);
        minimum->SetVariable(30,"m13",lines[1], step[4]);
        minimum->SetVariable(31,"q13",lines[1], step[0]);
        minimum->SetVariable(32,"m14",lines[1], step[4]);
        minimum->SetVariable(33,"q14",lines[1], step[0]);
        minimum->SetVariable(34,"m15",lines[1], step[4]);
        minimum->SetVariable(35,"q15",lines[1], step[0]);
        minimum->SetVariable(36,"m16",lines[1], step[4]);
        minimum->SetVariable(37,"q16",lines[1], step[0]);

        minimum->SetVariableLimits(22,-0.03,0.03);
        minimum->SetVariableLimits(23,-5.,5.);
        minimum->SetVariableLimits(24,-0.03,0.03);
        minimum->SetVariableLimits(25,-5.,5.);
        minimum->SetVariableLimits(26,-0.03,0.03);
        minimum->SetVariableLimits(27,-5.,5.);
        minimum->SetVariableLimits(28,-0.03,0.03);
        minimum->SetVariableLimits(29,-5.,5.);
        minimum->SetVariableLimits(30,-0.03,0.03);
        minimum->SetVariableLimits(31,-5.,5.);
        minimum->SetVariableLimits(32,-0.03,0.03);
        minimum->SetVariableLimits(33,-5.,5.);
        minimum->SetVariableLimits(34,-0.03,0.03);
        minimum->SetVariableLimits(35,-5.,5.);
        minimum->SetVariableLimits(36,-0.03,0.03);
        minimum->SetVariableLimits(37,-5.,5.);





    }
 
    // do the minimization
    minimum->Minimize();
 
    const double *xs = minimum->X();
    

    
    //Save the residuals

    std::vector<double> residuals={xs[0]-x[1],xs[1]-x[2],xs[2]-x[3],xs[3]-dy[1],xs[4]-dy[2],xs[5]-dy[3] };
    std::cerr<<contatore_tracce<<"  "<<xs[0]<<"    "<<xs[1]<<"    "<<xs[2]<<"  iteration    "<<0<<std::endl;
    res_iteration.at(0)=xs[0];
    res_iteration.at(1)=xs[1];
    res_iteration.at(2)=xs[2];

    res_iteration.at(3)=xs[3];
    res_iteration.at(4)=xs[4];
    res_iteration.at(5)=xs[5];

    // expected minimum is 0

    return residuals;
    minimum->Clear();
}
////Identical part but for iteration=1
else {
        std::cerr<<"ciao"<<std::endl;
    

    //ROOT::Math::Minimizer* minimum =
    //   ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
 
    ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "");

    // set tolerance , etc...
    minimum->SetMaxFunctionCalls(100000000); // for Minuit/Minuit2
    minimum->SetMaxIterations(10000);  // for GSL
    minimum->SetTolerance(0.00000001);
    minimum->SetPrintLevel(1);
    minimum->SetStrategy(2);
  
 
    // create function wrapper for minimizer
    // a IMultiGenFunction type
    ROOT::Math::Functor f(&chi_square,6+2*n);
    double step[19] = {0.0001, 0.0001, 0.0001,0.00001, 0.0001, 0.001,0.0001, 0.0001, 0.0001,0.0001, 0.0001, 0.0001,0.0001, 0.0001, 0.0001,0.0001, 0.0001, 0.0001,0.0001};
    // starting point
 
    //double variable[3] = { 39.4, 282.8, 326.8};

    double lines[16] = { 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,};

    double value=0;
 
    minimum->SetFunction(f);
 
    // Set the free variables to be minimized !
    //x positions
   
    minimum->SetVariable(0,"x2",result_iter_0[0], step[0]);
   
    minimum->SetVariable(1,"x3",result_iter_0[1], step[1]);
    
    minimum->SetVariable(2,"x4",result_iter_0[2], step[2]);

    minimum->SetVariableLimits(0,result_iter_0[0]-2.,result_iter_0[0]+2.);
    minimum->SetVariableLimits(1,result_iter_0[1]-2.,result_iter_0[1]+2.);
    minimum->SetVariableLimits(2,result_iter_0[2]-2.,result_iter_0[2]+2.);

  

    //dy of trackers

   
   
    minimum->SetVariable(3,"dy2",result_iter_0[3], step[0]);
   
    minimum->SetVariable(4,"dy3",result_iter_0[4], step[1]);
   
    minimum->SetVariable(5,"dy4",result_iter_0[5], step[2]);

     minimum->SetVariableLimits(3,-0.9,+1.1);
    minimum->SetVariableLimits(4,-1.3,0.7);
    minimum->SetVariableLimits(5,-0.1,+1.9);

    
    //lines

    minimum->SetVariable(6,"m1",lines[1], step[4]);
    minimum->SetVariable(7,"q1",lines[1], step[4]);

    minimum->SetVariableLimits(6,-0.03,0.03);
    minimum->SetVariableLimits(7,-5.,5.);
   
    if(n>1) {
        minimum->SetVariable(8,"m2",lines[1], step[4]);
        minimum->SetVariable(9,"q2",lines[1], step[4]);

        minimum->SetVariableLimits(8,-0.03,0.03);
        minimum->SetVariableLimits(9,-5.,5.);
    }
    if(n>2) {
        minimum->SetVariable(10,"m3",lines[1], step[4]);
        minimum->SetVariable(11,"q3",lines[1], step[4]);
        minimum->SetVariable(12,"m4",lines[1], step[4]);
        minimum->SetVariable(13,"q4",lines[1], step[4]);

        minimum->SetVariableLimits(10,-0.03,0.03);
        minimum->SetVariableLimits(11,-5.,5.);
        minimum->SetVariableLimits(12,-0.03,0.03);
        minimum->SetVariableLimits(13,-5.,5.);
    }
    if(n>4) {
        minimum->SetVariable(14,"m5",lines[1], step[4]);
        minimum->SetVariable(15,"q5",lines[1], step[4]);
        minimum->SetVariable(16,"m6",lines[1], step[4]);
        minimum->SetVariable(17,"q6",lines[1], step[4]);
        minimum->SetVariable(18,"m7",lines[1], step[4]);
        minimum->SetVariable(19,"q7",lines[1], step[4]);
        minimum->SetVariable(20,"m8",lines[1], step[4]);
        minimum->SetVariable(21,"q8",lines[1], step[4]);

        minimum->SetVariableLimits(14,-0.03,0.03);
        minimum->SetVariableLimits(15,-5.,5.);
        minimum->SetVariableLimits(16,-0.03,0.03);
        minimum->SetVariableLimits(17,-5.,5.);
        minimum->SetVariableLimits(18,-0.03,0.03);
        minimum->SetVariableLimits(19,-5.,5.);
        minimum->SetVariableLimits(20,-0.03,0.03);
        minimum->SetVariableLimits(21,-5.,5.);
    }
    if(n>8) {
        minimum->SetVariable(22,"m9",lines[1], step[4]);
        minimum->SetVariable(23,"q9",lines[1], step[4]);
        minimum->SetVariable(24,"m10",lines[1], step[4]);
        minimum->SetVariable(25,"q10",lines[1], step[4]);
        minimum->SetVariable(26,"m11",lines[1], step[4]);
        minimum->SetVariable(27,"q11",lines[1], step[4]);
        minimum->SetVariable(28,"m12",lines[1], step[4]);
        minimum->SetVariable(29,"q12",lines[1], step[4]);
        minimum->SetVariable(30,"m13",lines[1], step[4]);
        minimum->SetVariable(31,"q13",lines[1], step[4]);
        minimum->SetVariable(32,"m14",lines[1], step[4]);
        minimum->SetVariable(33,"q14",lines[1], step[4]);
        minimum->SetVariable(34,"m15",lines[1], step[4]);
        minimum->SetVariable(35,"q15",lines[1], step[4]);
        minimum->SetVariable(36,"m16",lines[1], step[4]);
        minimum->SetVariable(37,"q16",lines[1], step[4]);

        minimum->SetVariableLimits(22,-0.03,0.03);
        minimum->SetVariableLimits(23,-5.,5.);
        minimum->SetVariableLimits(24,-0.03,0.03);
        minimum->SetVariableLimits(25,-5.,5.);
        minimum->SetVariableLimits(26,-0.03,0.03);
        minimum->SetVariableLimits(27,-5.,5.);
        minimum->SetVariableLimits(28,-0.03,0.03);
        minimum->SetVariableLimits(29,-5.,5.);
        minimum->SetVariableLimits(30,-0.03,0.03);
        minimum->SetVariableLimits(31,-5.,5.);
        minimum->SetVariableLimits(32,-0.03,0.03);
        minimum->SetVariableLimits(33,-5.,5.);
        minimum->SetVariableLimits(34,-0.03,0.03);
        minimum->SetVariableLimits(35,-5.,5.);
        minimum->SetVariableLimits(36,-0.03,0.03);
        minimum->SetVariableLimits(37,-5.,5.);





    }
 
    // do the minimization
    minimum->Minimize();
 
    const double *xs = minimum->X();
    //std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "," << xs[2] << "): "
             // << minimum->MinValue()  << std::endl;

    

    //std::vector<double> residuals={x[1]-xs[0],x[2]-xs[1],x[3]-xs[2]};
    std::vector<double> residuals={xs[0]-x[1],xs[1]-x[2],xs[2]-x[3],xs[3]-dy[1],xs[4]-dy[2],xs[5]-dy[3] };
    std::cerr<<contatore_tracce<<"  "<<xs[0]<<"    "<<xs[1]<<"    "<<xs[2]<<"  iteration    "<<1<<std::endl;
 
    // expected minimum is 0

    
    return residuals;

}

 }



int main() {

    //it0 refers to iteration 0 and it1 for iteration 1

     std::ifstream i_lines("gaussian_output_lines_y_dy.txt");
     std::ifstream i_data("gaussian_output_hit_y_dy.txt");

     std::ofstream o1_it0(TString (namefile+"x1_it0_residuals.txt"));
     std::ofstream o2_it0(TString (namefile+"x2_it0_residuals.txt"));
     std::ofstream o3_it0(TString (namefile+"x3_it0_residuals.txt"));

     std::ofstream o1_it1(TString (namefile+"x1_it1_residuals.txt"));
     std::ofstream o2_it1(TString (namefile+"x2_it1_residuals.txt"));
     std::ofstream o3_it1(TString (namefile+"x3_it1_residuals.txt"));

     std::ofstream o_y1_it0(TString (namefile+"y1_it0_residuals.txt"));
     std::ofstream o_y2_it0(TString (namefile+"y2_it0_residuals.txt"));
     std::ofstream o_y3_it0(TString (namefile+"y3_it0_residuals.txt"));

     std::ofstream o_y1_it1(TString (namefile+"y1_it1_residuals.txt"));
     std::ofstream o_y2_it1(TString (namefile+"y2_it1_residuals.txt"));
     std::ofstream o_y3_it1(TString (namefile+"y3_it1_residuals.txt"));
    /////////////////////Initialization for iteration 0
    Double_t a=0;
    int cont=0;
    Double_t m_all[8000], q_all[8000], y_all[32000];
    while(i_lines>>a) {
        
        //std::cerr<<a<<std::endl;
        q_all[cont]=a;
        i_lines>>a;
        
        m_all[cont]=a;
        cont++;
    }

    cont=0;
    while(i_data>>a) {
    i_data>>a;
    y_all[cont]=a;
    cont++;
    }

    //for(int k=0;k<4000*n;k=k+4*n) {
    for(int k=0;k<2400*n;k=k+4*n) {  //->loop for all data. Max 4000 for n=8 and 2000 for n=16
        for(int l=0;l<n;l++) {
        dato a;
        a.m=m_all[(k+4*l)/4];
        a.q=q_all[(k+4*l)/4];
        a.y[0]=y_all[k+4*l];
        a.y[1]=y_all[k+4*l+1];
        a.y[2]=y_all[k+4*l+2];
        a.y[3]=y_all[k+4*l+3];
        dataset.push_back(a);
        }
        
        
        //cout of residuals
        
        std::vector<double> resid=NumericalMinimization();
        std::cerr<< dataset.size()<<std::endl;
        dataset.clear();
    
        o1_it0<<resid.at(0)<<std::endl;
        x2_it0.push_back(resid.at(0)+x[1]);
        o2_it0<<resid.at(1)<<std::endl;
        x3_it0.push_back(resid.at(1)+x[2]);
        o3_it0<<resid.at(2)<<std::endl;
        x4_it0.push_back(resid.at(2)+x[3]);

        o_y1_it0<<resid.at(3)<<std::endl;
        y2_it0.push_back(resid.at(3)+dy[1]);
        o_y2_it0<<resid.at(4)<<std::endl;
        y3_it0.push_back(resid.at(4)+dy[2]);
        o_y3_it0<<resid.at(5)<<std::endl;
        y4_it0.push_back(resid.at(5)+dy[3]);
    }
    o1_it0.close();
    o2_it0.close();
    o3_it0.close();
    o_y1_it0.close();
    o_y2_it0.close();
    o_y3_it0.close();
    i_data.close();
    i_lines.close();


    ///////////////////////////////////////////Iteration 1
    result_iter_0[0]=mean_v(x2_it0);
    result_iter_0[1]=mean_v(x3_it0);
    result_iter_0[2]=mean_v(x4_it0);
    result_iter_0[3]=mean_v(y2_it0);
    result_iter_0[4]=mean_v(y3_it0);
    result_iter_0[5]=mean_v(y4_it0);
    std::cerr<<result_iter_0[0]<<"  "<<result_iter_0[1]<<"  "<<result_iter_0[2]<<"  "<<result_iter_0[3]<<"  "<<result_iter_0[4]<<"  "<<result_iter_0[5]<<std::endl;



    std::ifstream i_lines2("gaussian_output_lines_y_dy.txt");
    std::ifstream i_data2("gaussian_output_hit_y_dy.txt");


    cont=0;

    while(i_lines2>>a) {
        //std::cerr<<a<<std::endl;
        q_all[cont]=a;
        i_lines2>>a;
        
        m_all[cont]=a;
        cont++;
    }

    cont=0;
    while(i_data2>>a) {
    i_data2>>a;
    y_all[cont]=a;
    cont++;
    }

    //for(int k=0;k<4000*n;k=k+4*n) {
    for(int k=0;k<2400*n;k=k+4*n) {
        for(int l=0;l<n;l++) {
        dato a;
        a.m=m_all[(k+4*l)/4];
        a.q=q_all[(k+4*l)/4];
        a.y[0]=y_all[k+4*l];
        a.y[1]=y_all[k+4*l+1];
        a.y[2]=y_all[k+4*l+2];
        a.y[3]=y_all[k+4*l+3];
        dataset.push_back(a);
        }
        
        
        //cout of residuals
        
        std::vector<double> resid=NumericalMinimization(1);
        std::cerr<< dataset.size()<<std::endl;
        dataset.clear();
    
        o1_it1<<resid.at(0)<<std::endl;
       
        o2_it1<<resid.at(1)<<std::endl;
        
        o3_it1<<resid.at(2)<<std::endl;
      

        o_y1_it1<<resid.at(3)<<std::endl;
       
        o_y2_it1<<resid.at(4)<<std::endl;
       
        o_y3_it1<<resid.at(5)<<std::endl;
        
    }





    return 0;


}


///////////////////////Plots


TH2F * draw_histo_2d_it0(int cont) {
    double a=0, b=0;
    std::ifstream i1(namefile+"x1_it0_residuals.txt");
    std::ifstream i2(namefile+"x2_it0_residuals.txt");
    std::ifstream i3(namefile+"x3_it0_residuals.txt");
    std::ifstream i_y1(namefile+"y1_it0_residuals.txt");
    std::ifstream i_y2(namefile+"y2_it0_residuals.txt");
    std::ifstream i_y3(namefile+"y3_it0_residuals.txt");

    if(cont==1) {
    TH2F * x_y_1=new TH2F("","y2 vs y3",200,-1.,1.,200,-1.,1.);
    while(i_y1>>a) {
        i_y2>>b;
        x_y_1->Fill(a,b);
    }
    i_y1.close();
    i_y2.close();
    return x_y_1;
    }
    if(cont==2) {
    TH2F * x_y_2=new TH2F("","y2 vs y4",200,-1.,1.,200,-1.,1.);
    while(i_y1>>a) {
        i_y3>>b;
        x_y_2->Fill(a,b);
    }
    i_y1.close();
    i_y3.close();
    return x_y_2;
    }
    if(cont==3) {
    TH2F * x_y_3=new TH2F("","y3 vs y4",200,-1.,1.,200,-1.,1.);
    while(i_y2>>a) {
        i_y3>>b;
        x_y_3->Fill(a,b);
    }
    i3.close();
    i_y3.close();
    return x_y_3;
    }
}

TH2F * draw_histo_2d_it1(int cont) {
    double a=0, b=0;
    std::ifstream i1(namefile+"x1_it1_residuals.txt");
    std::ifstream i2(namefile+"x2_it1_residuals.txt");
    std::ifstream i3(namefile+"x3_it1_residuals.txt");
    std::ifstream i_y1(namefile+"y1_it1_residuals.txt");
    std::ifstream i_y2(namefile+"y2_it1_residuals.txt");
    std::ifstream i_y3(namefile+"y3_it1_residuals.txt");

    if(cont==1) {
    TH2F * x_y_1=new TH2F("","x2 vs y2",200,-1.,1.,200,-1.,1.);
    while(i1>>a) {
        i_y1>>b;
        x_y_1->Fill(a,b);
    }
    i1.close();
    i_y1.close();
    return x_y_1;
    }
    if(cont==2) {
    TH2F * x_y_2=new TH2F("","x3 vs y3",200,-1.,1.,200,-1.,1.);
    while(i2>>a) {
        i_y2>>b;
        x_y_2->Fill(a,b);
    }
    i2.close();
    i_y2.close();
    return x_y_2;
    }
    if(cont==3) {
    TH2F * x_y_3=new TH2F("","x4 vs y4",200,-1.,1.,200,-1.,1.);
    while(i3>>a) {
        i_y3>>b;
        x_y_3->Fill(a,b);
    }
    i3.close();
    i_y3.close();
    return x_y_3;
    }
}


TH1F * draw_histo_it0(int cont) {
    double a=0;
    std::ifstream i1(namefile+"x1_it0_residuals.txt");
    std::ifstream i2(namefile+"x2_it0_residuals.txt");
    std::ifstream i3(namefile+"x3_it0_residuals.txt");
    std::ifstream i_y1(namefile+"y1_it0_residuals.txt");
    std::ifstream i_y2(namefile+"y2_it0_residuals.txt");
    std::ifstream i_y3(namefile+"y3_it0_residuals.txt");
    //x plot
    if(cont==1) {
    TH1F * x1=new TH1F("","Residuals x2",200,-1.,1.);
    while(i1>>a)
        x1->Fill(a);
    i1.close();
    return x1;
    }
    if(cont==2) {
    TH1F * x2=new TH1F("","Residuals x3",200,-1.,1.);
    while(i2>>a)
        x2->Fill(a);
    i2.close();
    return x2;
    }
    if(cont==3) {
    TH1F * x3=new TH1F("","Residuals x4",200,-1.,1.);
    while(i3>>a)
        x3->Fill(a);
    i3.close();
    return x3;
    }

    //y plot
     if(cont==4) {
    TH1F * y1=new TH1F("","Residuals y2",200,-1.,1.);
    while(i_y1>>a)
        y1->Fill(a);
    i_y1.close();
    return y1;
    }
    if(cont==5) {
    TH1F * y2=new TH1F("","Residuals y3",200,-1.,1.);
    while(i_y2>>a)
        y2->Fill(a);
    i_y2.close();
    return y2;
    }
    if(cont==6) {
    TH1F * y3=new TH1F("","Residuals y4",200,-1.,1.);
    while(i_y3>>a)
        y3->Fill(a);
    i_y3.close();
    return y3;
    }
    
    
}

TH1F * draw_histo_it1(int cont) {
    double a=0;
    std::ifstream i1(namefile+"x1_it1_residuals.txt");
    std::ifstream i2(namefile+"x2_it1_residuals.txt");
    std::ifstream i3(namefile+"x3_it1_residuals.txt");
    std::ifstream i_y1(namefile+"y1_it1_residuals.txt");
    std::ifstream i_y2(namefile+"y2_it1_residuals.txt");
    std::ifstream i_y3(namefile+"y3_it1_residuals.txt");
    //x plot
    if(cont==1) {
    TH1F * x1=new TH1F("","Residuals x2",200,-1.,1.);
    while(i1>>a)
        x1->Fill(a);
    i1.close();
    return x1;
    }
    if(cont==2) {
    TH1F * x2=new TH1F("","Residuals x3",200,-1.,1.);
    while(i2>>a)
        x2->Fill(a);
    i2.close();
    return x2;
    }
    if(cont==3) {
    TH1F * x3=new TH1F("","Residuals x4",200,-1.,1.);
    while(i3>>a)
        x3->Fill(a);
    i3.close();
    return x3;
    }

    //y plot
     if(cont==4) {
    TH1F * y1=new TH1F("","Residuals y2",200,-1.,1.);
    while(i_y1>>a)
        y1->Fill(a);
    i_y1.close();
    return y1;
    }
    if(cont==5) {
    TH1F * y2=new TH1F("","Residuals y3",200,-1.,1.);
    while(i_y2>>a)
        y2->Fill(a);
    i_y2.close();
    return y2;
    }
    if(cont==6) {
    TH1F * y3=new TH1F("","Residuals y4",200,-1.,1.);
    while(i_y3>>a)
        y3->Fill(a);
    i_y3.close();
    return y3;
    }
    
    
}


double dev_std(int cont) {
    std::vector<double> a;
    std::ifstream i1(namefile+"x1_residuals.txt");
    std::ifstream i2(namefile+"x2_residuals.txt");
    std::ifstream i3(namefile+"x3_residuals.txt");
    double b=0, dev=0, sum=0, mean=0;
    if(cont==1) {
        while(i1>>b) {
        a.push_back(b);
        //std::cerr<<b<<std::endl;
        sum+=b;
        }
    }
    if(cont==2) {
        while(i2>>b) {
        a.push_back(b);
        sum+=b;
        }
    }
    if(cont==3) {
        while(i3>>b) {
        a.push_back(b);
        sum+=b;
        }
    }
    
    mean=sum/a.size();

    sum=0;

    for(int i=0;i<a.size();i++) {
        sum+=pow(a.at(i)-mean,2);
    }
    
    dev=sqrt(sum/(a.size()-1));
    return dev;
}


void print_all_it0() {
  TH1F * x2=draw_histo_it0(1);
  TH1F * x3= draw_histo_it0(2);
  TH1F * x4 = draw_histo_it0(3);
  TH1F * y2=draw_histo_it0(4);
  TH1F * y3= draw_histo_it0(5);
  TH1F * y4 = draw_histo_it0(6);
  TH2F * x_y_2=draw_histo_2d_it0(1);
  TH2F * x_y_3 = draw_histo_2d_it0(2);
  TH2F * x_y_4 = draw_histo_2d_it0(3);
  std::cerr<<dev_std(1)<<"  "<<dev_std(2)<<"    "<<dev_std(3)<<std::endl;
  TF1 * g =new TF1("Gaussian","gaus",-0.5,0.5);
  /* x2->Fit(g);
  std::cerr<<g->GetParameter(2)<<std::endl;
  x3->Fit(g);
  std::cerr<<g->GetParameter(2)<<std::endl;
  x4->Fit(g);
  std::cerr<<g->GetParameter(2)<<std::endl;
  */

  TCanvas * c1 = new TCanvas("c1","c1",1600,1200);
  c1->Divide(3,3);
  c1->cd(1);
  //x2->GetXaxis()->SetRangeUser(-0.02,0.02);
  x2->GetXaxis()->SetTitle("x2 (cm)");
  x2->GetYaxis()->SetTitle("counts/(0.01 cm)");
  x2->Draw();
  c1->cd(2);
  //x3->GetXaxis()->SetRangeUser(-0.02,0.02);
  x3->GetXaxis()->SetTitle("x3 (cm)");
  x3->GetYaxis()->SetTitle("counts/(0.01 cm)");
  x3->Draw();
  c1->cd(3);
  //x4->GetXaxis()->SetRangeUser(0.02,0.02);
  x4->GetYaxis()->SetTitle("counts/(0.01 cm)");
  x4->GetXaxis()->SetTitle("x4 (cm)");
  x4->Draw();

  c1->cd(4);
  
  y2->GetXaxis()->SetTitle("y2 (cm)");
  y2->GetYaxis()->SetTitle("counts/(0.01 cm)");
  y2->Draw();
  c1->cd(5);
 
  y3->GetXaxis()->SetTitle("y3 (cm)");
  y3->GetYaxis()->SetTitle("counts/(0.01 cm)");
  y3->Draw();
  c1->cd(6);
  
  y4->GetYaxis()->SetTitle("counts/(0.01 cm)");
  y4->GetXaxis()->SetTitle("y4 (cm)");
  y4->Draw();

  c1->cd(7);
  
  x_y_2->GetXaxis()->SetTitle("x2 (cm)");
  x_y_2->GetYaxis()->SetTitle("y2 (cm)");
  x_y_2->Draw("COLZ");
  c1->cd(8);
 
  x_y_3->GetXaxis()->SetTitle("x3 (cm)");
  x_y_3->GetYaxis()->SetTitle("y3 (cm)");
  x_y_3->Draw("COLZ");
  c1->cd(9);
  
  x_y_4->GetYaxis()->SetTitle("y4 (cm)");
  x_y_4->GetXaxis()->SetTitle("x4 (cm)");
  x_y_4->Draw("COLZ");
  
}


void print_all_it1() {
  TH1F * x2=draw_histo_it1(1);
  TH1F * x3= draw_histo_it1(2);
  TH1F * x4 = draw_histo_it1(3);
  TH1F * y2=draw_histo_it1(4);
  TH1F * y3= draw_histo_it1(5);
  TH1F * y4 = draw_histo_it1(6);
  TH2F * x_y_2=draw_histo_2d_it1(1);
  TH2F * x_y_3 = draw_histo_2d_it1(2);
  TH2F * x_y_4 = draw_histo_2d_it1(3);
  std::cerr<<dev_std(1)<<"  "<<dev_std(2)<<"    "<<dev_std(3)<<std::endl;
  TF1 * g =new TF1("Gaussian","gaus",-0.5,0.5);
  /* x2->Fit(g);
  std::cerr<<g->GetParameter(2)<<std::endl;
  x3->Fit(g);
  std::cerr<<g->GetParameter(2)<<std::endl;
  x4->Fit(g);
  std::cerr<<g->GetParameter(2)<<std::endl;
  */

  TCanvas * c1 = new TCanvas("c1","c1",1600,1200);
  c1->Divide(3,3);
  c1->cd(1);
  //x2->GetXaxis()->SetRangeUser(-0.02,0.02);
  x2->GetXaxis()->SetTitle("x2 (cm)");
  x2->GetYaxis()->SetTitle("counts/(0.01 cm)");
  x2->Draw();
  c1->cd(2);
  //x3->GetXaxis()->SetRangeUser(-0.02,0.02);
  x3->GetXaxis()->SetTitle("x3 (cm)");
  x3->GetYaxis()->SetTitle("counts/(0.01 cm)");
  x3->Draw();
  c1->cd(3);
  //x4->GetXaxis()->SetRangeUser(0.02,0.02);
  x4->GetYaxis()->SetTitle("counts/(0.01 cm)");
  x4->GetXaxis()->SetTitle("x4 (cm)");
  x4->Draw();

  c1->cd(4);
  
  y2->GetXaxis()->SetTitle("y2 (cm)");
  y2->GetYaxis()->SetTitle("counts/(0.01 cm)");
  y2->Draw();
  c1->cd(5);
 
  y3->GetXaxis()->SetTitle("y3 (cm)");
  y3->GetYaxis()->SetTitle("counts/(0.01 cm)");
  y3->Draw();
  c1->cd(6);
  
  y4->GetYaxis()->SetTitle("counts/(0.01 cm)");
  y4->GetXaxis()->SetTitle("y4 (cm)");
  y4->Draw();

  c1->cd(7);
  
  x_y_2->GetXaxis()->SetTitle("x2 (cm)");
  x_y_2->GetYaxis()->SetTitle("y2 (cm)");
  x_y_2->Draw("COLZ");
  c1->cd(8);
 
  x_y_3->GetXaxis()->SetTitle("x3 (cm)");
  x_y_3->GetYaxis()->SetTitle("y3 (cm)");
  x_y_3->Draw("COLZ");
  c1->cd(9);
  
  x_y_4->GetYaxis()->SetTitle("y4 (cm)");
  x_y_4->GetXaxis()->SetTitle("x4 (cm)");
  x_y_4->Draw("COLZ");
  
}


void Scan() {

    TH2D * scan_2d=new TH2D("scan_2d","scan_2d",1000,-1.3,0.7,1000,-0.1,1.9);
    std::ifstream i_lines("gaussian_output_lines_y_dy.txt");
     std::ifstream i_data("gaussian_output_hit_y_dy.txt");
     std::ifstream i_lines_fit(TString (namefile+"lines_it0_fit.txt"));


    Double_t a=0;
    int cont=0;
    Double_t m_all[8000], q_all[8000], y_all[32000];
    while(i_lines_fit>>a) {
        
        //std::cerr<<a<<std::endl;
        q_all[cont]=a;
        i_lines_fit>>a;
        
        m_all[cont]=a;
        cont++;
    }

    cont=0;
    while(i_data>>a) {
    i_data>>a;
    y_all[cont]=a;
    cont++;
    }

    //for(int k=0;k<4000*n;k=k+4*n) {
    for(int k=0;k<4*n;k=k+4*n) {
        for(int l=0;l<n;l++) {
        dato a;
        a.m=m_all[(k+4*l)/4];
        a.q=q_all[(k+4*l)/4];
        a.y[0]=y_all[k+4*l];
        a.y[1]=y_all[k+4*l+1];
        a.y[2]=y_all[k+4*l+2];
        a.y[3]=y_all[k+4*l+3];
        dataset.push_back(a);
        }
    }



    for(int i=1;i<1001;i++) {
    for(int j=1;j<1001;j++) {
        
    double dy3=(scan_2d->GetXaxis())->GetBinCenter(i);
    double dy4=(scan_2d->GetYaxis())->GetBinCenter(j);
    //std::cerr<<dy3<<"   "<<dy4<<std::endl;

    double delta=0;
    for(int i=0; i<dataset.size();i++) {
       dato a=dataset.at(i);
       y[0]=a.y[0];
       y[1]=a.y[1];
       y[2]=a.y[2];
       y[3]=a.y[3];
       
       m=a.m;
       q=a.q;
       
       double sum1=0, sum2=0, sum3=0, sum4=0;
      

        sum1=(y[0]-q)*(y[0]-q)/(1*1);
        sum2=w2/w_tot * (y[1]-q-m*x[1]-0.1)*(y[1]-q-m*x[1]-0.1);
        sum3=w3/w_tot*(y[2]-q-m*x[2]+dy3)*(y[2]-q-m*x[2]+dy3);
        sum4=w4/w_tot*(y[3]-q-m*x[3]+dy4)*(y[3]-q-m*x[3]+dy4); 

         //std::cerr<<sum1<<" "<<sum2<<" "<<sum3<<" "<<sum4<<std::endl;

      

        
        //std::cerr<<log(sum1)<<" "<<log(sum2)<<" "<<log(sum3)<<" "<<log(sum4)<<std::endl;

        //std::vector<Double_t> sum_fin={sum1,sum2,sum3,sum4};
        //std::sort (sum_fin.begin(), sum_fin.end());

       
        delta=delta+sum1+sum2+sum3+sum4;
       
        

   }

    
    
    scan_2d->Fill(dy3,dy4,delta);
    }
    }
    scan_2d->GetZaxis()->SetRangeUser(0.,10.);
    scan_2d->Draw("COLZ");


}



