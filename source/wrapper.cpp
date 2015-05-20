/**
 * Minimal 2D Jacobi example. Code which is commented out demos how to
 * add a PPMWriter for output.
 */
/*
#include <libgeodecomp/io/simpleinitializer.h>
#include <libgeodecomp/io/ppmwriter.h>
#include <libgeodecomp/io/simplecellplotter.h>
#include <libgeodecomp/io/tracingwriter.h>
#include <libgeodecomp/parallelization/serialsimulator.h>
#include <libgeodecomp/parallelization/hiparsimulator.h>
*/
#include "../libgeodecomp/src/libgeodecomp.h"
#include <libflatarray/short_vec.hpp>
#include <cmath>


#define STENCIL_NONE 20
#define STENCIL_2D5PT 21 //  ! Default
#define STENCIL_2D9PT 22
#define STENCIL_3D7PT 23
#define STENCIL_3D27PT 24


/*
class MiniGhostContainer
{
public:
    static int numberOfVars;
};

int MiniGhostContainer::numberOfVars = 0;
*/
using namespace LibGeoDecomp;


class Cell
{
public:
    class API :
        public APITraits::HasFixedCoordsOnlyUpdate,
        public APITraits::HasUpdateLineX,
        //public APITraits::HasStencil<Stencils::VonNeumann<3, 1> >, 
	public APITraits::HasStencil<Stencils::Moore<3, 1> >, //contains all spatial neighbours
        public APITraits::HasOpaqueMPIDataType<Cell>,
        public APITraits::HasTorusTopology<3>,	// periodische randbedingung - cube ist konstant
        //public APITraits::HasCubeTopology<3>,
        //public APITraits::HasPredefinedMPIDataType<double>,
        public APITraits::HasSoA
    {};

    //	*** 2D5PT ***
    template<typename HOOD_OLD, typename HOOD_NEW>
    static void updateSingleStencil(HOOD_OLD& hoodOld, HOOD_NEW& hoodNew, int currentVar) {}

    template<typename HOOD_OLD, typename HOOD_NEW>
    static void updateLineX(HOOD_OLD& hoodOld, int indexEnd,
                            HOOD_NEW& hoodNew, int /* nanoStep */)
    {
        for ( ; hoodOld.index() < indexEnd; ++hoodOld.index(), ++hoodNew.index) {
            for(int currentVar = 0; currentVar < numberOfVars; ++currentVar) {
                updateSingleStencil(hoodOld, hoodNew, currentVar);
            }
        }

    }
    
    static int numberOfVars;  
    double temp[40] = {0};
};


LIBFLATARRAY_REGISTER_SOA(
    Cell,
    ((double)(temp)(40))
                          )

class Cell2D5PT
{
public:
    class API :
        public APITraits::HasFixedCoordsOnlyUpdate,
        public APITraits::HasUpdateLineX,
        //public APITraits::HasStencil<Stencils::VonNeumann<3, 1> >, 
	public APITraits::HasStencil<Stencils::Moore<3, 1> >, //contains all spatial neighbours
        public APITraits::HasOpaqueMPIDataType<Cell>,
        public APITraits::HasTorusTopology<3>,	// periodische randbedingung - cube ist konstant
        //public APITraits::HasCubeTopology<3>,
        //public APITraits::HasPredefinedMPIDataType<double>,
        public APITraits::HasSoA
    {};

    //	*** 2D5PT ***
    template<typename HOOD_OLD, typename HOOD_NEW>
    static void updateSingleStencil(HOOD_OLD& hoodOld, HOOD_NEW& hoodNew, int currentVar)
    {
        hoodNew.temp()[currentVar] =
            (hoodOld[FixedCoord<-1,  0,  0>()].temp()[currentVar] +
            hoodOld[FixedCoord< 0, -1,  0>()].temp()[currentVar] +
            hoodOld[FixedCoord< 0,  0,  0>()].temp()[currentVar] +
            hoodOld[FixedCoord< 0,  1,  0>()].temp()[currentVar] +
            hoodOld[FixedCoord< 1,  0,  0>()].temp()[currentVar]) * (1.0 / 5.0);
    }


    template<typename HOOD_OLD, typename HOOD_NEW>
    static void updateLineX(HOOD_OLD& hoodOld, int indexEnd,
                            HOOD_NEW& hoodNew, int /* nanoStep */)
    {
        for ( ; hoodOld.index() < indexEnd; ++hoodOld.index(), ++hoodNew.index) {
            for(int currentVar = 0; currentVar < numberOfVars; ++currentVar) {
                updateSingleStencil(hoodOld, hoodNew, currentVar);
            }
        }

    }
    
    static int numberOfVars;  
    double temp[40] = {0};
};


LIBFLATARRAY_REGISTER_SOA(
    Cell2D5PT,
    ((double)(temp)(40))
                          )


class Cell2D9PT
{
public:
    class API :
        public APITraits::HasFixedCoordsOnlyUpdate,
        public APITraits::HasUpdateLineX,
        //public APITraits::HasStencil<Stencils::VonNeumann<3, 1> >, 
	public APITraits::HasStencil<Stencils::Moore<3, 1> >, //contains all spatial neighbours
        public APITraits::HasOpaqueMPIDataType<Cell>,
        public APITraits::HasTorusTopology<3>,	// periodische randbedingung - cube ist konstant
        //public APITraits::HasCubeTopology<3>,
        //public APITraits::HasPredefinedMPIDataType<double>,
        public APITraits::HasSoA
    {};

    // 	*** 2D9PT ***
    template<typename HOOD_OLD, typename HOOD_NEW>
    static void updateSingleStencil(HOOD_OLD& hoodOld, HOOD_NEW& hoodNew, int currentVar)
    {
        hoodNew.temp()[currentVar] =
            (hoodOld[FixedCoord<-1, -1,  0>()].temp()[currentVar] +
            hoodOld[FixedCoord<-1,  0,  0>()].temp()[currentVar] +
            hoodOld[FixedCoord<-1,  1,  0>()].temp()[currentVar] +
		        
            hoodOld[FixedCoord< 0, -1,  0>()].temp()[currentVar] +
            hoodOld[FixedCoord< 0,  0,  0>()].temp()[currentVar] +
            hoodOld[FixedCoord< 0,  1,  0>()].temp()[currentVar] +
			
            hoodOld[FixedCoord< 1, -1,  0>()].temp()[currentVar] +
            hoodOld[FixedCoord< 1,  0,  0>()].temp()[currentVar] +
            hoodOld[FixedCoord< 1,  1,  0>()].temp()[currentVar]) * (1.0 / 9.0);			
    }
    
    template<typename HOOD_OLD, typename HOOD_NEW>
    static void updateLineX(HOOD_OLD& hoodOld, int indexEnd,
                            HOOD_NEW& hoodNew, int /* nanoStep */)
    {
        for ( ; hoodOld.index() < indexEnd; ++hoodOld.index(), ++hoodNew.index) {
            for(int currentVar = 0; currentVar < numberOfVars; ++currentVar) {
                updateSingleStencil(hoodOld, hoodNew, currentVar);
            }
        }

    }
    
    static int numberOfVars;  
    double temp[40] = {0};
};


LIBFLATARRAY_REGISTER_SOA(
    Cell2D9PT,
    ((double)(temp)(40))
                          )

class Cell3D7PT
{
public:
    class API :
        public APITraits::HasFixedCoordsOnlyUpdate,
        public APITraits::HasUpdateLineX,
        //public APITraits::HasStencil<Stencils::VonNeumann<3, 1> >, 
	public APITraits::HasStencil<Stencils::Moore<3, 1> >, //contains all spatial neighbours
        public APITraits::HasOpaqueMPIDataType<Cell>,
        public APITraits::HasTorusTopology<3>,	// periodische randbedingung - cube ist konstant
        //public APITraits::HasCubeTopology<3>,
        //public APITraits::HasPredefinedMPIDataType<double>,
        public APITraits::HasSoA
    {};
    
    //  *** 3D7PT ***
    template<typename HOOD_OLD, typename HOOD_NEW>
    static void updateSingleStencil(HOOD_OLD& hoodOld, HOOD_NEW& hoodNew, int currentVar)
    {
        hoodNew.temp()[currentVar] =
            (hoodOld[FixedCoord<-1,  0,  0>()].temp()[currentVar] +		
            hoodOld[FixedCoord< 0, -1,  0>()].temp()[currentVar] +
            hoodOld[FixedCoord< 0,  0,  0>()].temp()[currentVar] +
            hoodOld[FixedCoord< 0,  1,  0>()].temp()[currentVar] +	
            hoodOld[FixedCoord< 1,  0,  0>()].temp()[currentVar] +
            hoodOld[FixedCoord< 0,  0, -1>()].temp()[currentVar] +
            hoodOld[FixedCoord< 0,  0,  1>()].temp()[currentVar]) / (1.0 / 7.0); 
    }

    template<typename HOOD_OLD, typename HOOD_NEW>
    static void updateLineX(HOOD_OLD& hoodOld, int indexEnd,
                            HOOD_NEW& hoodNew, int /* nanoStep */)
    {
        for ( ; hoodOld.index() < indexEnd; ++hoodOld.index(), ++hoodNew.index) {
            for(int currentVar = 0; currentVar < numberOfVars; ++currentVar) {
                updateSingleStencil(hoodOld, hoodNew, currentVar);
            }
        }

    }
    
    static int numberOfVars;  
    double temp[40] = {0};
};


LIBFLATARRAY_REGISTER_SOA(
    Cell3D7PT,
    ((double)(temp)(40))
                          )


class Cell3D27PT
{
public:
    class API :
        public APITraits::HasFixedCoordsOnlyUpdate,
        public APITraits::HasUpdateLineX,
        //public APITraits::HasStencil<Stencils::VonNeumann<3, 1> >, 
	public APITraits::HasStencil<Stencils::Moore<3, 1> >, //contains all spatial neighbours
        public APITraits::HasOpaqueMPIDataType<Cell>,
        public APITraits::HasTorusTopology<3>,	// periodische randbedingung - cube ist konstant
        //public APITraits::HasCubeTopology<3>,
        //public APITraits::HasPredefinedMPIDataType<double>,
        public APITraits::HasSoA
    {};
   
    //  *** 3D27PT ***
    template<typename HOOD_OLD, typename HOOD_NEW>
    static void updateSingleStencil(HOOD_OLD& hoodOld, HOOD_NEW& hoodNew, int currentVar)
    {
        hoodNew.temp()[currentVar] = 
            // *** BACK ***
            (hoodOld[FixedCoord<-1, -1, -1>()].temp()[currentVar] +
            hoodOld[FixedCoord<-1,  0, -1>()].temp()[currentVar] +
            hoodOld[FixedCoord<-1,  1, -1>()].temp()[currentVar] +
            
            hoodOld[FixedCoord< 0, -1, -1>()].temp()[currentVar] +
            hoodOld[FixedCoord< 0,  0, -1>()].temp()[currentVar] +
            hoodOld[FixedCoord< 0,  1, -1>()].temp()[currentVar] +
            
            hoodOld[FixedCoord< 1, -1, -1>()].temp()[currentVar] +
            hoodOld[FixedCoord< 1,  0, -1>()].temp()[currentVar] +
            hoodOld[FixedCoord< 1,  1, -1>()].temp()[currentVar] +
            
            // *** MIDDLE ***
            hoodOld[FixedCoord<-1, -1,  0>()].temp()[currentVar] +
            hoodOld[FixedCoord<-1,  0,  0>()].temp()[currentVar] +
            hoodOld[FixedCoord<-1,  1,  0>()].temp()[currentVar] +
            
            hoodOld[FixedCoord< 0, -1,  0>()].temp()[currentVar] +
            hoodOld[FixedCoord< 0,  0,  0>()].temp()[currentVar] +
            hoodOld[FixedCoord< 0,  1,  0>()].temp()[currentVar] +
            
            hoodOld[FixedCoord< 1, -1,  0>()].temp()[currentVar] +
            hoodOld[FixedCoord< 1,  0,  0>()].temp()[currentVar] +
            hoodOld[FixedCoord< 1,  1,  0>()].temp()[currentVar] +
            
            // *** FRONT ***
            hoodOld[FixedCoord<-1, -1,  1>()].temp()[currentVar] +
            hoodOld[FixedCoord<-1,  0,  1>()].temp()[currentVar] +
            hoodOld[FixedCoord<-1,  1,  1>()].temp()[currentVar] +
            
            hoodOld[FixedCoord< 0, -1,  1>()].temp()[currentVar] +
            hoodOld[FixedCoord< 0,  0,  1>()].temp()[currentVar] +
            hoodOld[FixedCoord< 0,  1,  1>()].temp()[currentVar] +
            
            hoodOld[FixedCoord< 1, -1,  1>()].temp()[currentVar] +
            hoodOld[FixedCoord< 1,  0,  1>()].temp()[currentVar] +
            hoodOld[FixedCoord< 1,  1,  1>()].temp()[currentVar]) / (1.0 / 27.0);
    }

    template<typename HOOD_OLD, typename HOOD_NEW>
    static void updateLineX(HOOD_OLD& hoodOld, int indexEnd,
                            HOOD_NEW& hoodNew, int /* nanoStep */)
    {
        for ( ; hoodOld.index() < indexEnd; ++hoodOld.index(), ++hoodNew.index) {
            for(int currentVar = 0; currentVar < numberOfVars; ++currentVar) {
                updateSingleStencil(hoodOld, hoodNew, currentVar);
            }
        }

    }
    
    static int numberOfVars;  
    double temp[40] = {0};
};


LIBFLATARRAY_REGISTER_SOA(
    Cell3D27PT,
    ((double)(temp)(40))
                          )


template<typename CELL>
class CellInitializer : public SimpleInitializer<CELL>
{
public:
    using SimpleInitializer<CELL>::Topology;

    CellInitializer(const unsigned dimX, const unsigned dimY, const unsigned dimZ, const unsigned num_timesteps, int numVars, double *sourceTotal_, double *spikes_ , int *spikeLoc) :
    SimpleInitializer<CELL>(Coord<3>(dimX, dimY, dimZ), num_timesteps),
    dimX(dimX),
    dimY(dimY),
    dimZ(dimZ),
    numberOfVars(numVars)
    {
        sourceTotal = sourceTotal_;
        spikes = spikes_;
        spikeLocation = spikeLoc;
    }

    virtual void grid(GridBase<CELL, 3> *ret)
    {
        
        //Cell::numberOfVars = numberOfVars;
        		
        CoordBox<3> rect = ret->boundingBox();        
        for (CoordBox<3>::Iterator i = rect.begin(); i != rect.end(); ++i)
        {
//coord c;
//Topology.normalize(c)
//seed = c.toIndex(dim)

//neue frage: wenn wir nach je koordinate eine pseudorand erzeugen, dann sind alle pro variable gleich also vlt sowas wie
//seed = global_gridsize * currentVar + c.toIndex(dim)
//oder wie ist das mit dem seed und der randomfolge??


/*Tabelle fuer: welcher rang welche start und end coords hat (global)
fuer ein paar problemgroessen
um checkerboard zu testen
*/
//new function in town - pseudoRand(*i);
	    //seede den random nr generator 
	    
            CELL cell = CELL();

            for( int currentVar = 0; currentVar < numberOfVars; currentVar++ )
            {
                cell.temp[currentVar] = 1;//Random::gen_d();
            }
            
            ret->set(*i, cell);
	}
/*	
	// set multiple first spikes
        Coord<3> c( (unsigned) spikeLocation[1],
	            (unsigned) spikeLocation[2],
		    (unsigned) spikeLocation[3]);
		
        if (rect.inBounds(c)) 
        {
	    CELL cell = CELL();
 
            for( int currentVar = 0; currentVar < numberOfVars; currentVar++ )
            {
                cell.temp[currentVar] = spikes[currentVar];

                std::cout << "WARNING---------------------------------------------------\n"
	              << "WARNING: We're at Location " <<  spikeLocation[1] << ", " << spikeLocation[2] << ", " << spikeLocation[3] << "\n"
	              << "WARNING: and we 're setting the initial spike " << std::setprecision (15) << spikes[currentVar] << " into the grid\n"
		      << "WARNING---------------------------------------------------\n";
            }
            
            ret->set(c, cell);

        }
*/		    
        // update sourceTotal on every node
	for( int currentVar = 0; currentVar < numberOfVars; currentVar++ )
        {
	    sourceTotal[currentVar] = sourceTotal[currentVar] + spikes[currentVar];
        }

    }

private:
    unsigned dimX;
    unsigned dimY;
    unsigned dimZ;
    int numberOfVars;
    double *sourceTotal;
    double *spikes;
    int *spikeLocation;
};

template<typename CELL>
class ToleranceChecker : public Clonable<ParallelWriter<CELL>, ToleranceChecker<CELL>>
{
public:
    using typename ParallelWriter<CELL>::RegionType;
    using typename ParallelWriter<CELL>::CoordType;
    typedef typename ParallelWriter<CELL>::GridType GridType;
	
    ToleranceChecker(const unsigned outputPeriod = 1, int numberOfVars = 1, double errorTol = 2500.0,  double *sourceTotal = NULL) :
    Clonable<ParallelWriter<CELL>, ToleranceChecker>("", outputPeriod),
    num_vars(numberOfVars),
    err_tol(errorTol)
    {
        src_total = sourceTotal;
    }
	
    void stepFinished(
		const GridType& grid,
		const RegionType& validRegion,
		const CoordType& globalDimensions,
		unsigned step,
		WriterEvent event,
		std::size_t rank,
		bool lastCall)
    {		
	// summing the local grid (not sure what access pattern is better)
	for( int currentVar = 0 ; currentVar < num_vars ; ++currentVar )
	{
	    for (typename RegionType::Iterator i = validRegion.begin(); i != validRegion.end(); ++i) { 
	        localSum[currentVar] += grid.get(*i).temp[currentVar];
	    }
	} 
			
	// so we're done summing the local grid
        if(lastCall)
        {
            //tag fuer einzigartigkeit
            MPI_Allreduce(localSum, globalSum, num_vars, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                        
            // we're still missing the feature how many of our grids (variables in this case) should be summed
            //(see MG_BUFINIT.F)
            
            // checking error tolerance
            for( int currentVar = 0 ; currentVar < num_vars ; ++currentVar )
            {
                if( rank == 0 )
            	    std::cout << "globalSum(" << step << ") = " << std::setprecision (15) << globalSum[currentVar] << "\n";

            	if( ( std::abs(src_total[currentVar] - globalSum[currentVar]) / src_total[currentVar] ) > err_tol )
            	    std::cout << "error_tol(" << step << ") has not been met for Variable " << currentVar << "\n";
                    
                // reset local and global sum since we're done with this step (better memcpy)
                localSum[currentVar] = 0;
                globalSum[currentVar] = 0;
            }
            
        }
    }
        
private: 
    double localSum[40] = {0};
    double globalSum[40] = {0}; 
    int num_vars;
    double err_tol;
    double *src_total;
};

template<typename CELL>
class SpikeJab : public Steerer<CELL>
{
public:
    using typename Steerer<CELL>::CoordType;
    using typename Steerer<CELL>::GridType;
    using typename Steerer<CELL>::Topology;

    SpikeJab (const unsigned ioPeriod, int numVars = 1, int numSpikes = 1, double *sourceTotal_ = NULL, double *spikes_ = NULL, int *spikeLocation_ = NULL) :
        Steerer<CELL>(ioPeriod),
        numberOfVars(numVars),
        numberOfSpikes(numSpikes)
    {
	sourceTotal = sourceTotal_;
	spikes = spikes_;
	spikeLocation = spikeLocation_;
    }

    void nextStep(
        GridType *grid,
        const Region<Topology::DIM>& validRegion,
        const CoordType& globalDimensions,
        unsigned step,
        SteererEvent event,
        std::size_t rank,
        bool lastCall,
        typename Steerer<CELL>::SteererFeedback *feedback)
    {
        // gibt es einen cooleren weg dem steerer zu sagen dass er den allerletzten step nichtmehr machen braucht?
        if ( step == numberOfSpikes * this->getPeriod() ) 
	{
	    return;
	}
		
	currentSpike = ( (int)step / this->getPeriod() );
	
	const Coord<3> currentCoord( spikeLocation[(currentSpike * 4) + 1], 
				     spikeLocation[(currentSpike * 4) + 2],
				     spikeLocation[(currentSpike * 4) + 3]);
	
	// setting spikes even in ghostzones							
	if(validRegion.count(currentCoord))
	{
            //not really necessary is it?
	    //Cell cell = grid->get(currentCoord);

            CELL cell = CELL();
	    for( int currentVar = 0; currentVar < numberOfVars; currentVar++ )
            {
                cell.temp[currentVar] = spikes[(currentSpike * numberOfVars) + currentVar];  

                std::cout << "WARNING---------------------------------------------------\n"
                          << "WARNING: We're at Location " <<  spikeLocation[(currentSpike * 4) + 1] << ", " << spikeLocation[(currentSpike * 4) + 2] << ", " << spikeLocation[(currentSpike * 4) + 3] << "\n"
                          << "WARNING: We're rank: " << rank << "\n"
                          << "WARNING: and we 're setting spike " << std::setprecision (15) << cell.temp[currentVar] << " into grid at timestep " << step << "\n"
                          << "WARNING---------------------------------------------------\n";
            }
            grid->set(currentCoord, cell);	
			
	}
	// everybody refresh sourceTotal - whether or not it belongs to this local grid
	if( lastCall )
	{
	    for( int currentVar = 0; currentVar < numberOfVars; currentVar++ )
	    {
	    	sourceTotal[currentVar] = sourceTotal[currentVar] + spikes[(currentSpike * numberOfVars) + currentVar];
	    }
	}
		  
    }

private:
    int numberOfVars;
    int numberOfSpikes;
    double *sourceTotal;
    double *spikes;
    int *spikeLocation;
    int currentSpike;
};

int Cell::numberOfVars = 0;

//LINE DOMINANT VS ROW DOMINANT - WATCH OUT !
//double spikes[*num_spikes][*num_vars], double spike_loc[*num_spikes][4]
extern "C" void simulate_(int *nx, int *ny, int *nz, int *stencil, int *num_vars, int *num_spikes, int *num_tsteps, double *err_tol, double *source_total, double *spikes, int *spike_loc)
{
    Cell::numberOfVars = *num_vars;

    std::cout << " <<<<<<<< " << *stencil << " >>>>>>>>>\n";
        
    switch (*stencil) {
        case STENCIL_NONE  : // do somethin
                             break;
        case STENCIL_2D5PT : // do somethin else
                             break;
        case STENCIL_2D9PT : // do somethin else but almost like the others
                             break;
        case STENCIL_3D7PT : // freakin' do it 
                             break;
        case STENCIL_3D27PT: // i swear u need to do it here too
                             break;
        default : std::cerr << " not a stecil \n";  break;
    }
 
//    Typemaps::initializeMaps();
    {
	//SerialSimulator<Cell> sim(new CellInitializer(dimX, dimY, num_timesteps));
	CellInitializer<Cell> *init =  new CellInitializer<Cell>(*nx, *ny, *nz, (*num_tsteps * *num_spikes), *num_vars, source_total, spikes, spike_loc);
	
	//CheckerBoarding		
	HiParSimulator::HiParSimulator<Cell, RecursiveBisectionPartition<3> > *sim = 
                        new HiParSimulator::HiParSimulator<Cell, RecursiveBisectionPartition<3> >(
			init,
			//MPILayer().rank() ? 0 : new TracingBalancer(new NoOpBalancer()),
			0,
			*num_tsteps,
			1);
	
	ToleranceChecker<Cell> *toleranceChecker = new ToleranceChecker<Cell>(1, *num_vars, *err_tol, source_total);
	sim->addWriter(toleranceChecker);

	SpikeJab<Cell> *spikeJab = new SpikeJab<Cell>(*num_tsteps, *num_vars, *num_spikes, source_total, spikes, spike_loc);
	sim->addSteerer(spikeJab);
	
	sim->run();
    }
    
}
