
#include<iostream>
#include<fstream>
#include "TRandom3.h"
#include "TPolyLine3D.h"
#include "TView3D.h"
#include "TGraph2D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TPad.h"
#include "TStyle.h"
#include "TMinuit.h"

using namespace std;
std::ofstream o1("gaussian_output_hit_y_dy.txt");
std::ofstream o2("gaussian_output_lines_y_dy.txt");
std::ofstream o3("gaussian_output_hit_z_dz.txt");
std::ofstream o4("gaussian_output_lines_z_dz.txt");

//GEOMETRIC SYSTEM OF REFERENCE
//x->beam direction
//y->vertical dimension
//z-> direction perpendicular to the beam and parallel to the ground

//Origin is set in the middle of T1

//Measurements in cm
TRandom3 *r=new TRandom3(42);



struct hit{
    double x=0, y=0, z=0; 
};

hit SetHit(hit a) {
    hit b;
    b.x=a.x;
    b.y=a.y;
    b.z=a.z;
    return b;
}





//Class for a single tracker
class tracker {
    public:
    double GetD_X() {return d_x;};
    double GetD_Y() {return d_y;};
    double GetD_Z() {return d_z;};
    void SetD_X(double a) {d_x=a;};
    void SetD_Y(double a) {d_y=a;};
    void SetD_Z(double a) {d_z=a;};
    void Fill_Histo(hit a) {
                            hit_yz->Fill(a.y-d_y,a.z-d_z);
                            };
    void Fill_Hit(hit a) {
                            hit_vector_yz.push_back(a);
                        };
    //True case
    void Fill_Histo_True(hit a) {
                            hit_yz_true->Fill(a.y-d_y,a.z-d_z);
                            };
    void Fill_Hit_True(hit a) {
                            hit_vector_yz_true.push_back(a);
                        };

    void SetD(double a, double b, double c) {d_x=a;
                                            d_y=b;
                                            d_z=c;};
    double GetL_X() {return l_x;};
    double GetL_Y() {return l_y;};
    TH2F * GetHisto() {return hit_yz;};
    vector<hit> GetHits() {return hit_vector_yz;};
    TH2F * GetHistoTrue() {return hit_yz_true;};
    vector<hit> GetHitsTrue() {return hit_vector_yz_true;};
    std::vector<hit> GetVertex() {return vertex;};
    void SetL_X(double a) {l_x=a;};
    void SetL_Y(double a) {l_y=a;};
    void SetVertex() {
                    hit a;
                    a.x=d_x; a.y=d_y-5; a.z=d_z+5;
                    vertex.push_back(a);
                    a.x=d_x; a.y=d_y-5; a.z=d_z-5;
                    vertex.push_back(a);
                    a.x=d_x; a.y=d_y+5; a.z=d_z-5;
                    vertex.push_back(a);
                    a.x=d_x; a.y=d_y+5; a.z=d_z+5;
                    vertex.push_back(a);
                    };
    

    protected:
    double d_x=0, d_y=0, d_z=0; //position in GEOMETRIC setup
    double l_x=10., l_y=10.;  //dimension in GEOMETRIC system of reference---->different from ITS OWN
    TH2F * hit_yz = new TH2F("hit_yz","",312,-5,5,312,-5,5); // map of hits centered on tracker center
    TH2F * hit_yz_true = new TH2F("hit_yz_true","",312,-5,5,312,-5,5); // map of hits centered on tracker center (true point)
    vector<hit> vertex;
    vector<hit> hit_vector_yz;
    vector<hit> hit_vector_yz_true;
};

//Class for a line in 3d
class line {
    public:
    double GetPhi() {return phi;};
    double GetTheta() {return theta;};
    hit GetP0() {return p0;};
    hit Projection(tracker tr) {
                                hit track_hit;
                                track_hit.x=tr.GetD_X();
                                track_hit.y=p0.y+track_hit.x*tan(theta)-tr.GetD_Y();
                                track_hit.z=p0.z+track_hit.x*tan(phi)/sin(theta)-tr.GetD_Z();                
                                return track_hit;
                                 };
    void SetPhi(double a) {phi=a;};
    void SetTheta(double a) {theta=a;};
    void Setx0(double a) {p0.x=a;};
    void Sety0(double a) {p0.y=a;};
    void Setz0(double a) {p0.z=a;};
    void LineInit();
    bool IsGood(tracker tr);

    protected:
    double theta=0, phi=0;
    hit p0;
};

void line::LineInit() {
    double a=r->Uniform(-5.,5.);
    p0.y=a;
    a=r->Uniform(-5.,5.);
    p0.z=a;
    //a=r->Uniform(-M_PI/20,M_PI/20);
    a=r->Uniform(-0.07,0.07);
    theta=a;
    //a=r->Uniform(-0.01,0.01);
    a=r->Uniform(-0.07,0.07);
    phi=a;
}




//It contains the 4 tracker shifted randomly
class setup {
    public:

    void Initialize();
    vector<line> sim_data(int n);
    void draw_geom(line l1);
    tracker GetTracker(int a) {
                            if(a==1)
                                return t1;
                            if(a==2)
                                return t2;
                            if(a==3)
                                return t3;
                            if(a==4)
                                return t4;
                             };

    protected:
    tracker t1, t2, t3, t4;


};

void setup::Initialize() {
    double delta_x=0, delta_y=0, delta_z=0;
    //T2
    delta_x=39.4; //Measured 39.4
    delta_y=0.1; //Measured 0.1
    delta_z=0.4; //Unknown
    t2.SetD(delta_x,delta_y,delta_z);

    //T3
    delta_x=282.8; //Measured 282.8
    delta_y=-0.3; //Measured -0.3
    delta_z=-0.2; //Unknown
    t3.SetD(delta_x,delta_y,delta_z);

    //T4
    delta_x=326.8; //Measured 326.8
    delta_y=0.9; //Measured 0.9
    delta_z=-0.5; //Unknown
    t4.SetD(delta_x,delta_y,delta_z);
    t1.SetVertex();
    t2.SetVertex();
    t3.SetVertex();
    t4.SetVertex();
}


void setup::draw_geom(line l1) {
//Tracker 1
    TPolyLine3D *contour_1 = new TPolyLine3D(0);
    contour_1->SetPoint(0,t1.GetD_X(),t1.GetD_Y()-5.,t1.GetD_Z()+5.);
    contour_1->SetPoint(1,t1.GetD_X(),t1.GetD_Y()-5.,t1.GetD_Z()-5.);
    contour_1->SetPoint(2,t1.GetD_X(),t1.GetD_Y()+5.,t1.GetD_Z()-5.);
    contour_1->SetPoint(3,t1.GetD_X(),t1.GetD_Y()+5.,t1.GetD_Z()+5.);
    contour_1->SetPoint(4,t1.GetD_X(),t1.GetD_Y()-5.,t1.GetD_Z()+5.);

    //Tracker 2
    TPolyLine3D *contour_2 = new TPolyLine3D(0);
    contour_2->SetPoint(0,t2.GetD_X(),t2.GetD_Y()-5.,t2.GetD_Z()+5.);
    contour_2->SetPoint(1,t2.GetD_X(),t2.GetD_Y()-5.,t2.GetD_Z()-5.);
    contour_2->SetPoint(2,t2.GetD_X(),t2.GetD_Y()+5.,t2.GetD_Z()-5.);
    contour_2->SetPoint(3,t2.GetD_X(),t2.GetD_Y()+5.,t2.GetD_Z()+5.);
    contour_2->SetPoint(4,t2.GetD_X(),t2.GetD_Y()-5.,t2.GetD_Z()+5.);

    //Tracker 3
    TPolyLine3D *contour_3 = new TPolyLine3D(0);
    contour_3->SetPoint(0,t3.GetD_X(),t3.GetD_Y()-5.,t3.GetD_Z()+5.);
    contour_3->SetPoint(1,t3.GetD_X(),t3.GetD_Y()-5.,t3.GetD_Z()-5.);
    contour_3->SetPoint(2,t3.GetD_X(),t3.GetD_Y()+5.,t3.GetD_Z()-5.);
    contour_3->SetPoint(3,t3.GetD_X(),t3.GetD_Y()+5.,t3.GetD_Z()+5.);
    contour_3->SetPoint(4,t3.GetD_X(),t3.GetD_Y()-5.,t3.GetD_Z()+5.);

    //Tracker 4
    TPolyLine3D *contour_4 = new TPolyLine3D(0);
    contour_4->SetPoint(0,t4.GetD_X(),t4.GetD_Y()-5.,t4.GetD_Z()+5.);
    contour_4->SetPoint(1,t4.GetD_X(),t4.GetD_Y()-5.,t4.GetD_Z()-5.);
    contour_4->SetPoint(2,t4.GetD_X(),t4.GetD_Y()+5.,t4.GetD_Z()-5.);
    contour_4->SetPoint(3,t4.GetD_X(),t4.GetD_Y()+5.,t4.GetD_Z()+5.);
    contour_4->SetPoint(4,t4.GetD_X(),t4.GetD_Y()-5.,t4.GetD_Z()+5.);

    //Line
    TPolyLine3D *muon_line = new TPolyLine3D(0);
    for(double i=0;i<350;i++) {
        double a=0, b=0, c=0;
        a=l1.GetP0().x+i*cos(l1.GetTheta())*cos(l1.GetPhi());
        b=l1.GetP0().y+i*sin(l1.GetTheta())*cos(l1.GetPhi());
        c=l1.GetP0().z+i*sin(l1.GetPhi());
        //cerr<<l1.GetTheta()<<"  "<<l1.GetPhi()<<endl;
        //cerr<<a<<"  "<<b<<" "<<c<<std::endl;
        if(a>-10 && a<350 && abs(b)<25 && abs(c)<25)
            muon_line->SetPoint(i,a,b,c);
        else   
            break;
    }
    muon_line->SetLineColor(kRed);
    TH3F * geom_view=new TH3F("geom_view","geom_view",1,-5,350,1,-25,25,1,-20,20);
    gStyle->SetOptStat(0);
    geom_view->Draw();
    contour_1->Draw("SAME");
    contour_2->Draw("SAME");
    contour_3->Draw("SAME");
    contour_4->Draw("SAME");
    muon_line->Draw("SAME");
}


//Needed for IsGood function
bool FindPosition(hit a, tracker tr) {
    vector<hit> vertex2=tr.GetVertex(); 

    for(int i=0;i<vertex2.size();i++) {
        double pos=0;
        hit p, q;
        if(i==vertex2.size()-1) {
            p=SetHit(vertex2.at(i));
            q=SetHit(vertex2.at(0));
        }
        else {
            p=SetHit(vertex2.at(i));
            q=SetHit(vertex2.at(i+1));
        }
        pos=((a.z-p.z)*(q.y-p.y)-(a.y-p.y)*(q.z-p.z));
        //cerr<<pos<<endl;
        if(pos>0) continue;
        if(pos<0) return false;
    }
    return true;
    
}

//Find if the line passes through a polygon in which all vertices have the same x coordinate
bool line::IsGood(tracker tr) {
    hit track_hit;
    track_hit=Projection(tr);
    return FindPosition(track_hit,tr);

}
//VERY IMPORTANT: Change the true hit in a tracker with a one shifted randomly with a gaussian in y and z
hit tracker_resolution (hit true_point) {
    hit out_point;
    double a=r->Gaus(0,0.003);
    double b=r->Gaus(0,0.003);
    //double a=0, b=0;
    out_point.x=true_point.x;
    out_point.y=true_point.y+a;
    out_point.z=true_point.z+b;
    return out_point;
}




//Find n muons that pass through the 4 trackers
vector<line> setup::sim_data(int n) {

    vector<line> all_lines;
    int line_counter=0;
    

    while(line_counter<n) {
        //Define a line
        line l1;
    

        //Check if line passes through the 4 trackers
        bool rep=false;
        bool rep1=false, rep2=false, rep3=false, rep4=false;
        while(rep==false) {
            l1.LineInit(); //Random parameters

            hit track_hit;
            rep1=l1.IsGood(t1);
            //cerr<<"tracker1 "<<rep1<<endl;
            if(rep1==false)
                continue;
            rep2=l1.IsGood(t2);
            //cerr<<"tracker2 "<<rep2<<endl;
            if(rep2==false)
                continue;
            rep3=l1.IsGood(t3);
            //cerr<<"tracker3 "<<rep3<<endl;
            if(rep3==false)
                continue;
            rep4=l1.IsGood(t4);
            if(rep4==false)
                continue;
            //cerr<<line_counter<<endl;
            //cerr<<rep1<<"   "<<rep2<<"  "<<rep3<<"  "<<rep4<<endl;
            rep=rep1&&rep2&&rep3&&rep4;  
        }
        all_lines.push_back(l1);   
        line_counter++;
        if(line_counter%10==0)
            cerr<<line_counter<<endl; 

        //Fill with hits the trackers. They are corrected based upon dy and dz
        
        hit point_tracker_true=l1.Projection(t1);
        hit point_tracker=tracker_resolution(point_tracker_true);
        t1.Fill_Hit(point_tracker);
        //t1.Fill_Histo(point_tracker);
        t1.Fill_Hit_True(point_tracker_true);
        //t1.Fill_Histo_True(point_tracker_true);

        point_tracker_true=l1.Projection(t2);
        point_tracker=tracker_resolution(point_tracker_true);
        t2.Fill_Hit(point_tracker);
        //t2.Fill_Histo(point_tracker);
        t2.Fill_Hit_True(point_tracker_true);
        //t2.Fill_Histo_True(point_tracker_true);

        point_tracker_true=l1.Projection(t3);
        point_tracker=tracker_resolution(point_tracker_true);
        t3.Fill_Hit(point_tracker);
        //t3.Fill_Histo(point_tracker);
        t3.Fill_Hit_True(point_tracker_true);
        //t3.Fill_Histo_True(point_tracker_true);


        point_tracker_true=l1.Projection(t4);
        point_tracker=tracker_resolution(point_tracker_true);
        t4.Fill_Hit(point_tracker);
        //t4.Fill_Histo(point_tracker);
        t4.Fill_Hit_True(point_tracker_true);
        //t4.Fill_Histo_True(point_tracker_true);

        


        //if(line_counter==6)
        //    draw_geom(l1);    

    }
    return all_lines;

}



void simulation() {
    setup s1;
    //Define a random setup by changing the Delta around the measured value
    //Values are considered as a difference from first tracker
    s1.Initialize();
    //Generation of n muons that pass through the 4 trackers
    vector<line> lines=s1.sim_data(8000);
    cerr<<lines.size()<<endl;
    for(int i=0;i<lines.size();i++) {
        o2<<lines.at(i).GetP0().y<<"    "<<tan(lines.at(i).GetTheta())<<endl<<endl;
        o1<<s1.GetTracker(1).GetHits().at(i).x<<"   "<<s1.GetTracker(1).GetHits().at(i).y<<endl;
        o1<<s1.GetTracker(2).GetHits().at(i).x<<"   "<<s1.GetTracker(2).GetHits().at(i).y<<endl;
        o1<<s1.GetTracker(3).GetHits().at(i).x<<"   "<<s1.GetTracker(3).GetHits().at(i).y<<endl;
        o1<<s1.GetTracker(4).GetHits().at(i).x<<"   "<<s1.GetTracker(4).GetHits().at(i).y<<endl<<endl;

        o4<<lines.at(i).GetP0().z<<"    "<<tan(lines.at(i).GetPhi())/cos(lines.at(i).GetTheta())<<endl<<endl;
        o3<<s1.GetTracker(1).GetHits().at(i).x<<"   "<<s1.GetTracker(1).GetHits().at(i).z<<endl;
        o3<<s1.GetTracker(2).GetHits().at(i).x<<"   "<<s1.GetTracker(2).GetHits().at(i).z<<endl;
        o3<<s1.GetTracker(3).GetHits().at(i).x<<"   "<<s1.GetTracker(3).GetHits().at(i).z<<endl;
        o3<<s1.GetTracker(4).GetHits().at(i).x<<"   "<<s1.GetTracker(4).GetHits().at(i).z<<endl<<endl;
    }

    //tracker tr_ex=s1.GetTracker(1);
    //tr_ex.GetHisto()->Draw("COLZ");
    
    //return s1;
}