#include "utils.h"


const Eigen::Matrix3d v2t(vector<double>& vec){

    double pi = M_PI;
    double x = vec[0];
    double y = vec[1];
    double heading = vec[2];
    
    if(heading < 0){
        heading += 2*pi;
    }
    else if(heading > pi*2){
        heading -= pi*2;
    }

    Eigen::Matrix3d T(3, 3);
    T<<cos(heading), -sin(heading), x, sin(heading), cos(heading), y, 0, 0, 1;

    return T; 
}

const vector<double> t2v(const Eigen::Ref<Eigen::Matrix3d>& mtx){

    double pi = M_PI;

    double x = mtx(0,2);
    double y = mtx(1,2);
    double heading = atan2(mtx(1,0), mtx(0,0));

    if(heading < 0){
        heading += 2*pi;
    }
    else if(heading > pi*2){
        heading -= pi*2;
    }

    vector<double> vec = {x, y, heading};

    return vec;    
}

Gaussian multGaussians(const Gaussian& g1, const Gaussian& g2){


    //cout<<"Gaissans multiplication - Start"<<endl;

    // Null hypothesis case
    if(g1.isNull && g2.isNull){
        return Gaussian(true, g1.weight*g2.weight);
    }

    else if(!g1.isNull and g2.isNull){
        return g1;
    }

    else if(g1.isNull and !g2.isNull){
        return g1;
    }


    else{
        int dim = g1.dim;

        double pi = M_PI;

        // multiplication of two multivariate Gaussian - diagonal case
        vector<double> inv_cov1 = {1.0/g1.cov[0], 1.0/g1.cov[1], 1.0/g1.cov[2]};
        vector<double> inv_cov2 = {1.0/g2.cov[0], 1.0/g2.cov[1], 1.0/g2.cov[2]};

        // covariance of result Gaussian (diagonal only)
        vector<double> cov = {(g1.cov[0]*g2.cov[0])/(g1.cov[0]+g2.cov[0]),
                                   (g1.cov[1]*g2.cov[1])/(g1.cov[1]+g2.cov[1]),
                                   (g1.cov[2]*g2.cov[2])/(g1.cov[2]+g2.cov[2])};

        vector<double> mean1 = {g1.mean[0], g1.mean[1], g1.mean[2]};
        vector<double> mean2 = {g2.mean[0], g2.mean[1], g2.mean[2]};

        // angle correction
        if(mean1[2] < pi/2 && mean2[2] > 3.0/2.0*pi){
            mean2[2] -= pi*2;
        }

        else if(mean2[2] < pi/2 && mean1[2] > 3.0/2.0*pi){
            mean1[2] -= pi*2;
        }

        vector<double> mean = {cov[0]*inv_cov1[0]*mean1[0] + cov[0]*inv_cov2[0]*mean2[0],
                                    cov[1]*inv_cov1[1]*mean1[1] + cov[1]*inv_cov2[1]*mean2[1],
                                    cov[2]*inv_cov1[2]*mean1[2] + cov[2]*inv_cov2[2]*mean2[2]};

        return Gaussian(dim, g1.weight*g2.weight, mean, cov);

    }


    // cout<<"Gaissans multiplication - End"<<endl;

    // cout<<"############ DEBUG RESULT ############"<<endl;

    // cout<<g->mean[0]<<" "<<g->mean[1]<<" "<<g->mean[2]<<endl;
    // cout<<g->cov[0]<<" "<<g->cov[1]<<" "<<g->cov[2]<<endl;



    // cout<<"############ DEBUG RESULT ############"<<endl;

}

double calcNormalPDF(vector<double> x, vector<double> mean, vector<double> cov, int dim){
    
    double pi = M_PI;
    
    double det = (cov[0]*cov[1]*cov[2]);

    // calc denominator
    double denominator = pow(pow(2*pi, (double)dim)*det, 0.5);

    // calc mahalanobis distance
    double dx = x[0] - mean[0];
    double dy = x[1] - mean[1];

    double a = x[2];
    double a_mean = mean[2];

    //// angle correction : clip [0, 2*pi]
    if(a<0){
        a += 2*pi;
    }
    else if(a>2*pi){
        a -= 2*pi;
    }

    if(a_mean<0){
        a_mean += 2*pi;
    }
    else if(a_mean>2*pi){
        a_mean -= 2*pi;
    }

    if(a < pi/2 && a_mean > 3.0/2.0*pi){
            a_mean -= pi*2;
    }

    else if(a_mean < pi/2 && a > 3.0/2.0*pi){
        a -= pi*2;
    }

    double da = a - a_mean;
    //double da = x[2] - mean[2];




    double mahalanobis_dist = (dx*dx/cov[0] + dy*dy/cov[1] + da*da/cov[2]);
    
    // calc numerator
    double numerator = exp(-0.5*(mahalanobis_dist));
    double probability = numerator/denominator;    

    return probability;

}

void showMatrix(Eigen::Matrix3d _T){
    cout<<_T(0, 0)<<" "<<_T(0, 1)<<" "<<_T(0, 2)<<endl;
    cout<<_T(1, 0)<<" "<<_T(1, 1)<<" "<<_T(1, 2)<<endl;
    cout<<_T(2, 0)<<" "<<_T(2, 1)<<" "<<_T(2, 2)<<endl;
}

vector<vector<int>> cartProduct (const vector<vector<int>>& v) {
    vector<vector<int>> s = {{}};
    for (auto& u : v) {
        vector<vector<int>> r;
        for (auto& x : s) {
            for (auto y : u) {
                r.push_back(x);
                r.back().push_back(y);
            }
        }
        s.swap(r);
    }
    return s;
}


vector<Gaussian> exactSampling(vector<vector<Gaussian>>& mixtures, int dim, bool reparam){
    
    // make combinations - start
    if(mixtures.size()==1){
        return mixtures[0];
    }

    vector<vector<int>> sizes(mixtures.size());

    for(int i=0; i<mixtures.size(); i++){
        vector<int> size_temp(mixtures[i].size());
        for(int j=0; j<mixtures[i].size(); j++){
            size_temp[j] = j;
        }
        sizes[i] = move(size_temp);
    }

    vector<vector<int>> combs = cartProduct(sizes);

    // define parameters 

    int nb_gaussians = combs.size();

    int cnt = 0;

    vector<vector<double>> means(combs.size(), vector<double>(3));
    vector<vector<double>> covs(combs.size(), vector<double>(3));
    vector<double> ws(combs.size());


    // exact sampling
    //omp_set_num_threads(5);
    //#pragma omp parallel for
    //for(auto& comb : combs){
    for(int c=0; c<combs.size(); c++){
        Gaussian g_mul = move(multGaussians(mixtures[0][combs[c][0]], mixtures[1][combs[c][1]]));

        for(int i=2; i<mixtures.size(); i++){
            g_mul = (multGaussians(g_mul, mixtures[i][combs[c][i]]));

        }

        //g_mul.showInfo();
        double denominator = calcNormalPDF(g_mul.mean, g_mul.mean, g_mul.cov, dim)+1e-100;
        double numerator = 1.0;

        for(int i=0; i<mixtures.size(); i++){
            numerator *= calcNormalPDF(g_mul.mean, mixtures[i][combs[c][i]].mean, mixtures[i][combs[c][i]].cov, dim);
        }
        numerator += 1e-100;

        double w = (numerator/denominator);
        ws[cnt] = w;

        means[cnt] = (g_mul.mean);
        covs[cnt] = (g_mul.cov);

        cnt += 1;
    }

    

    double w_max = *max_element(ws.begin(), ws.end());
    
    vector<vector<double>> means_new;
    vector<vector<double>> covs_new;
    vector<double> ws_new;

    for(int k=0; k<means.size(); k++){
        if(ws[k] >= w_max*0.3){
            means_new.push_back(means[k]);
            covs_new.push_back(covs[k]);
            ws_new.push_back(ws[k]);
        }
    }

    means = move(means_new);
    covs = move(covs_new);
    ws = move(ws_new);


    if(reparam){
        //cout<<"Before : "<<means.size()<<endl;
        getMode(means, covs, ws, 10, 3, 300);
        //cout<<"After : "<<means.size()<<endl;
    }

    double w_sum = accumulate(ws.begin(), ws.end(), 0.0);

    for(auto& w : ws){
        w /= w_sum;
    }


    vector<vector<double>> mn = {means[0]};
    vector<vector<double>> cn = {covs[0]};
    vector<double> wn = {ws[0]};
    //////////////////////////////////////////////////////////////// EM here

    // vector<Gaussian> mixture_result(ws.size());

    // for(int i=0; i<ws.size(); i++){
    //     mixture_result[i] = (Gaussian(dim, ws[i], means[i], covs[i]));
    // }

    vector<Gaussian> mixture_result(wn.size());

    for(int i=0; i<wn.size(); i++){
        mixture_result[i] = (Gaussian(dim, wn[i], mn[i], cn[i]));
    }

    return mixture_result;
}

void normalizeVector(vector<double>& vec){

    double sum = accumulate(vec.begin(), vec.end(), 0.0);

    for(auto& v : vec){
        v /= sum;
    }
}

void printResults(vector<Point>& points, int num_points)
{
    int i = 0;
    printf("Number of points: %u\n"
        " x     y     z     cluster_id\n"
        "-----------------------------\n"
        , num_points);
    while (i < num_points)
    {
          printf("%5.2lf %5.2lf %5.2lf: %d\n",
                 points[i].x,
                 points[i].y, points[i].z,
                 points[i].clusterID);
          ++i;
    }
}

void getClass(vector<Point>& pts, vector<int>& _class){
    for(int i=0; i<pts.size(); i++){
        _class[i] = pts[i].clusterID;
    }
}


void getMode(vector<vector<double>>& means, vector<vector<double>>& covs, vector<double>& ws, int max_iter, int dim, int nb_bin){
    

    vector<vector<double>> xs = means; // copy

    double entropy_ref = 0.0;

    double minDelta=0.0;
    double maxDelta=0.0;
    double bin_size=0.0;

    // for(auto& w : ws){
    //     cout<<w<<" ";
    // }
    // cout<<endl;

    

    for(int iter=0; iter<max_iter; iter++){
        int n = means.size();

        vector<vector<double>> fx(n, vector<double>(3));
        
        for(int i=0; i<n; i++){
            
            // calculate ks(xi)
            vector<double> ks(n);

            for(int j=0; j<n; j++){
                ks[j] = ws[i]*calcNormalPDF(xs[j], means[i], covs[i], dim); // ks(xi) += k_ij
            }

            normalizeVector(ks);

            vector<double> kinfm_sum = {0, 0, 0};
            vector<double> kinfm_mean_sum = {0, 0, 0};

            // calculate f(xi)
            for(int j=0; j<n; j++){                
                for(int k=0; k<dim; k++){
                    kinfm_sum[k] += ks[j]*(1/(covs[j][k]+1e-100));
                    kinfm_mean_sum[k] += ks[j]*(1/(covs[j][k]+1e-100))*xs[j][k];
                }
            }

            for(int k=0; k<dim; k++){
                fx[i][k] = (1.0/(kinfm_sum[k]+1e-100))*(kinfm_mean_sum[k]) - xs[i][k];
            }
        }
        vector<double> delta(n);
        
        for(int i=0; i<n; i++){
            for(int k=0; k<3; k++){
                xs[i][k] += fx[i][k];
                delta[i] += sqrt(pow(fx[i][k], 2));
            }
        }
        // stop criteria
        if(iter == 0){
            minDelta = 0.0;
            maxDelta = *max_element(delta.begin(), delta.end());
            bin_size = (maxDelta - minDelta)/((double)(nb_bin-1));   
        }
        vector<double> hist(nb_bin, 0);

        for(auto& d : delta){
            double idx = (int)(trunc((d - minDelta)/bin_size));
            if(idx >= hist.size()){
                idx = hist.size()-1;
            }
            else if(idx < 0){
                idx = 0;
            }

            hist[(int)(idx)] += 1.0;
        }
        double entropy = 0.0;

        for(auto& h : hist){
            h /= ((double)n);
            if(h != 0){
                entropy += (-h*log(h));
            }
        }
        double entropy_def = abs(entropy - entropy_ref);

        if(entropy_def < 1e-7){

            //cout<<iter<<" th iter : Done "<<entropy_def<<endl;
            break;
        }
        entropy_ref = entropy;

    }

    int MINIMUM_POINTS = 1;
    double EPSILON = 1;

    vector<Point> pts;
    vec2pts(xs, pts);

    DBSCAN ds(MINIMUM_POINTS, EPSILON, pts);
    ds.run();

    vector<int> class_id(pts.size());

    getClass(ds.m_points, class_id);


    // print_std_vector(class_id);

    int nb_cluster = *max_element(class_id.begin(), class_id.end());

    for(auto& ii : class_id){
        if(ii == -1){
            cout<<"db scan is not working..."<<endl;
        }
    }

    vector<vector<double>> means_new(nb_cluster, vector<double>(dim));
    vector<vector<double>> covs_new(nb_cluster, vector<double>(dim));
    vector<double> ws_new(nb_cluster);
    
    for(int i=0; i<class_id.size(); i++){
        ws_new[class_id[i]-1] += ws[i];
    }

    for(int i=0; i<class_id.size(); i++){
        for(int j=0; j<dim; j++){
            means_new[class_id[i]-1][j] += (ws[i]/ws_new[class_id[i]-1]*means[i][j]);
        }
    }

    for(int i=0; i<class_id.size(); i++){
        for(int j=0; j<dim; j++){
            covs_new[class_id[i]-1][j] += (ws[i]/ws_new[class_id[i]-1] + pow((means[i][j] - means_new[class_id[i]-1][j]), 2.0));
        }
    }

    means = means_new;
    covs = covs_new;
    ws = ws_new;

    double w_max = *max_element(ws.begin(), ws.end());

    vector<vector<double>> means_new2;
    vector<vector<double>> covs_new2;
    vector<double> ws_new2;

    for(int k=0; k<means.size(); k++){
        if(ws[k] >= w_max*0.3){
            means_new2.push_back(means[k]);
            covs_new2.push_back(covs[k]);
            ws_new2.push_back(ws[k]);
        }
    }

    means = move(means_new2);
    covs = move(covs_new2);
    ws = move(ws_new2);
}



void vec2pts(vector<vector<double>>& pts_in, vector<Point>& pts_out){

    for(int i=0; i<pts_in.size(); i++){
        pts_out.push_back(move(Point(pts_in[i])));
    }
}


vector<int> intersection(vector<int> &v1, vector<int> &v2){
    vector<int> v3;


    set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(v3));
    return v3;
}

void copy_vector(igraph_vector_t *v, vector<int>& vec){
    for(int i=0; i<vec.size(); i++){
        vec[i] = move((long int)VECTOR(*v)[i]);
    }
}


void print_vector(igraph_vector_t *v, FILE *f) {
    long int i;
    for (i = 0; i < igraph_vector_size(v); i++) {
        fprintf(f, " %li", (long int) VECTOR(*v)[i]);
    }
    fprintf(f, "\n");
}

void print_std_vector(vector<int>& vec){
    for(int i=0; i<vec.size(); i++){
        cout<<vec[i]<<" ";
    }
    cout<<endl;
}


void readCSV(string path, vector<vector<int>>& _idx_v, vector<vector<double>>& _odom){
    string data = path;
    ifstream in(data.c_str());

    typedef boost::tokenizer<boost::escaped_list_separator<char>> Tokenizer;
    vector<string> vec;
    string line;

    vector<vector<int>> idx_from_to;
    vector<vector<double>> odometry;

    while (getline(in,line))
    {
        Tokenizer tok(line);
        vec.assign(tok.begin(),tok.end());

        int idx_from = stoi(vec[1]);
        int idx_to = stoi(vec[2]);
        double x = stod(vec[3]);
        double y = stod(vec[4]);
        double h = stod(vec[5]);

        vector<int> idxs = {idx_from, idx_to};
        vector<double> odom = {x, y, h};

        _idx_v.push_back(idxs);
        _odom.push_back(odom);        
    }

}

void readCSV_MH(string path, vector<vector<int>>& _idx_v, vector<vector<vector<double>>>& _odom, int last_idx){
    string data = path;
    ifstream in(data.c_str());

    typedef boost::tokenizer<boost::escaped_list_separator<char>> Tokenizer;
    vector<string> vec;
    string line;

    vector<vector<int>> idx_from_to;
    vector<vector<double>> odometry;

    int cnt = 0;
    while (getline(in,line))
    {
        Tokenizer tok(line);
        vec.assign(tok.begin(),tok.end());

        int idx_from = stoi(vec[1]);
        int idx_to = stoi(vec[3]);
        int nb_odom = stoi(vec[5]);

        if(idx_to > last_idx){
            break;
        }

        if(nb_odom == 1){
            double x = stod(vec[6]);
            double y = stod(vec[7]);
            double h = stod(vec[8]);

            vector<int> idxs = {idx_from, idx_to};
            vector<vector<double>> odom = {{x, y, h}};

            _idx_v.push_back(idxs);
            _odom.push_back(odom);
        }
        else{
            double x = stod(vec[6]);
            double y = stod(vec[7]);
            double h = stod(vec[8]);
            double x2 = stod(vec[9]);
            double y2 = stod(vec[10]);
            double h2 = stod(vec[11]);

            vector<int> idxs = {idx_from, idx_to};
            vector<vector<double>> odom = {{x, y, h},{x2, y2, h2}};

            _idx_v.push_back(idxs);
            _odom.push_back(odom);
        }        
    }

}

void readCSV_MH_GT(string path, vector<vector<double>>& _odom, int last_idx){
    string data = path;
    ifstream in(data.c_str());

    typedef boost::tokenizer<boost::escaped_list_separator<char>> Tokenizer;
    vector<string> vec;
    string line;

    vector<vector<int>> idx_from_to;
    vector<vector<double>> odometry;
    int cnt = 0;

    while (getline(in,line))
    {
        if(cnt > last_idx){
            break;
        }
        Tokenizer tok(line);
        vec.assign(tok.begin(),tok.end());

        double x = stod(vec[0]);
        double y = stod(vec[1]);
        double h = stod(vec[2]);

        vector<double> odom = {x, y, h};

        _odom.push_back(odom);
        cnt += 1;
    }

}

double error_per(vector<double>& est, vector<double>& gt){
    return sqrt(pow(est[0]-gt[0], 2.0) + pow(est[1]-gt[1], 2));
}


double calcDist(const vector<double>& p1, const vector<double>& p2){
    double h1 = p1[2];
    double h2 = p2[2];

    double pi = M_PI;

    if(h1<0){
        h1 += 2*pi;
    }
    else if(h1>2*pi){
        h1 -= 2*pi;
    }

    if(h2<0){
        h2 += 2*pi;
    }
    else if(h2>2*pi){
        h2 -= 2*pi;
    }

    if(h1 < pi/2 && h2 > 3.0/2.0*pi){
            h2 -= pi*2;
    }

    else if(h2 < pi/2 && h1 > 3.0/2.0*pi){
        h1 -= pi*2;
    }

    return sqrt(pow(p1[0] - p2[0],2)+pow(p1[1] - p2[1],2)+pow(p1[2] - p2[2],2));
}


