#include "ocam_functions.h"
#include <Eigen/Core>
#include <Eigen/Geometry>

void create_virtualview_undistortion_LUT( 
        cv::Mat &mapx, cv::Mat &mapy, 
        struct ocam_model &ocam_model, 
        float azimuth, float elevation,
        float zoom )
{
     int i, j;
     int width = mapx.cols; //New width
     int height = mapx.rows;//New height     
     float *data_mapx = mapx.ptr<float>();
     float *data_mapy = mapy.ptr<float>();
     float Nxc = height/2.0;
     float Nyc = width/2.0;
     float Nz  = -width/zoom;
     double M[3];
     double m[2];

     Eigen::Quaterniond rot =
         Eigen::AngleAxisd( azimuth, Eigen::Vector3d::UnitZ() )
         * Eigen::AngleAxisd( elevation, Eigen::Vector3d::UnitY() );
     
     for (i=0; i<height; i++)
         for (j=0; j<width; j++)
         {   
             M[0] = (i - Nxc);
             M[1] = (j - Nyc);
             M[2] = Nz;

             // apply rotation
             Eigen::Map<Eigen::Vector3d> p( M );
             p = rot * p;

             world2cam(m, M, &ocam_model);
             *( data_mapx + i*width+j ) = (float) m[1];
             *( data_mapy + i*width+j ) = (float) m[0];
         }
}

int main(int argc, char *argv[])
{   
    if( argc < 3 )
    {
        std::cout << "usage: omniview <videofile> <calibfile>" << std::endl;
        exit(0);
    }

    std::string videofile = argv[1];
    std::string calibfile = argv[2];

    struct ocam_model model; 
    get_ocam_model(&model, calibfile.c_str());
    
    cv::VideoCapture vc( videofile ); 
    if( !vc.isOpened() )
    {
        std::cout << "could not open video file " << videofile << std::endl;
        exit(0);
    }
    double fps = vc.get(CV_CAP_PROP_FPS);
    int delay = 1000 / fps;

    std::cout <<
        "fps: " << fps << std::endl;

    cv::Mat frame;
    bool running = true, init = true;
    float azimuth = 0, 
          elevation = 0, 
          zoom = 2.0;

    cv::Size size( 800, 800 );
    cv::Mat view( size, CV_8UC3 );
    cv::Mat 
        mapx( size, CV_32FC1 ),
        mapy( size, CV_32FC1 );

    float elevation_step = M_PI * 10 / 180;
    float azimuth_step = M_PI * 10 / 180;
    float zoom_step = 0.5;

    while( running )
    {
        char key = cv::waitKey( delay );
        switch( key )
        {
            case 'q': running = false; break;
            case 'w': elevation += elevation_step; break;
            case 's': elevation -= elevation_step; break;
            case 'a': azimuth -= azimuth_step; break;
            case 'd': azimuth += azimuth_step; break;
            case 'r': zoom += zoom_step; break;
            case 'f': zoom -= zoom_step; break;
        }

        if( key != -1 || init )
        {
            azimuth = std::max( -(float)M_PI, std::min( azimuth, (float)M_PI ) );
            elevation = std::max( -(float)(M_PI/2.0), std::min( elevation, (float)(M_PI/2.0)) );
            zoom = std::max( 0.1f, std::min( 20.0f, zoom ) );

            create_virtualview_undistortion_LUT( 
                    mapx, mapy, model,
                    azimuth, elevation + M_PI/2.0, zoom );

            std::cout << "azimuth: "
                << azimuth / M_PI * 180
                << " elevation: "
                << elevation / M_PI * 180
                << " zoom: "
                << zoom 
                << std::endl;

            init = false;
        }

        vc >> frame;
        if( frame.empty() )
            break;

        cv::remap( frame, view, mapx, mapy, cv::INTER_LINEAR, cv::BORDER_CONSTANT, cv::Scalar(0,0,0) );

        cv::imshow( "omniview", view );
    }
}
