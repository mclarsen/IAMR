
#include <cstdio>

#include <AMReX_CArena.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Amr.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_BLProfiler.H>

#include <AMReX_AmrLevel.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Conduit_Blueprint.H>
#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint.hpp>
#include <conduit/conduit_relay.hpp>

#include <ascent.hpp>

using namespace amrex;

int
main (int   argc,
      char* argv[])
{
    std::cout<<"Hello!\n";
    amrex::Initialize(argc,argv);

    BL_PROFILE_REGION_START("main()");
    BL_PROFILE_VAR("main()", pmain);

    const Real run_strt = ParallelDescriptor::second();

    int  max_step;
    int  num_steps;
    Real strt_time;
    Real stop_time;

    ParmParse pp;

    max_step  = -1;
    num_steps = -1;
    strt_time =  0.0;
    stop_time = -1.0;

    pp.query("max_step",  max_step);
    pp.query("num_steps", num_steps);
    pp.query("strt_time", strt_time);
    pp.query("stop_time", stop_time);

    if (strt_time < 0.0)
    {
        amrex::Abort("MUST SPECIFY a non-negative strt_time");
    }

    if (max_step < 0 && stop_time < 0)
    {
        amrex::Abort("Exiting because neither max_step nor stop_time is non-negative.");
    }

    Amr* amrptr = new Amr;

    amrptr->init(strt_time,stop_time);

    if (num_steps > 0)
    {
        if (max_step < 0)
        {
            max_step = num_steps + amrptr->levelSteps(0);
        }
        else
        {
            max_step = std::min(max_step, num_steps + amrptr->levelSteps(0));
        }

	amrex::Print() << "Using effective max_step = " << max_step << '\n';
    }
    //
    // If we set the regrid_on_restart flag and if we are *not* going to take
    // a time step then we want to go ahead and regrid here.
    //
    if (amrptr->RegridOnRestart())
    {
        if (    (amrptr->levelSteps(0) >= max_step ) ||
                ( (stop_time >= 0.0) &&
                  (amrptr->cumTime() >= stop_time)  )    )
        {
            //
            // Regrid only!
            //
            amrptr->RegridOnly(amrptr->cumTime());
        }
    }

    while ( amrptr->okToContinue()                            &&
           (amrptr->levelSteps(0) < max_step || max_step < 0) &&
           (amrptr->cumTime() < stop_time || stop_time < 0.0) )
    {
        amrptr->coarseTimeStep(stop_time);

        std::cout<<"0000000000000000000000000000000000000000000000000000000000000000000000\n";
        //const std::list<std::string>& plot_vars = amrptr->statePlotVars();
        Vector<const MultiFab*> mfs;
        Vector<Geometry> geoms;
        Vector<int> level_steps;
        Vector<std::string> var_names;
        Vector<IntVect> ref_ratios;
        const int level = 0;

        const DescriptorList& desc_lst = amrptr->getLevel(0).get_desc_lst();
        std::vector<std::pair<int,int> > plot_var_map;
        for (int typ = 0; typ < desc_lst.size(); typ++)
        {
          for (int comp = 0; comp < desc_lst[typ].nComp();comp++)
          {
            if (amrptr->isStatePlotVar(desc_lst[typ].name(comp)) &&
                desc_lst[typ].getType() == IndexType::TheCellType())
            {
                plot_var_map.push_back(std::pair<int,int>(typ,comp));
                var_names.push_back(desc_lst[typ].name(comp));
            }
          }
        }

        for(int level = 0; level <= amrptr->finestLevel(); ++level)
        {

          // for each level
          const int nGrow = 1;
          const int n_data_items = plot_var_map.size();

          MultiFab*  level_mf = new MultiFab(amrptr->getLevel(level).boxArray(),
                                             amrptr->getLevel(level).DistributionMap(),
                                             n_data_items,
                                             nGrow,
                                             MFInfo(),
                                             amrptr->getLevel(level).Factory());
          mfs.push_back(level_mf);
          MultiFab* this_dat = 0;
          //
          // Cull data from state variables --
          //
          for (int i = 0; i < static_cast<int>(plot_var_map.size()); i++)
          {
            int type  = plot_var_map[i].first;
            int comp = plot_var_map[i].second;
            this_dat = &amrptr->getLevel(level).get_new_data(type);
            MultiFab::Copy(*level_mf,*this_dat,comp,i,1,nGrow);
          }

          geoms.push_back(amrptr->getLevel(level).Geom());
          level_steps.push_back(amrptr->levelSteps(level));
          ref_ratios.push_back(amrptr->getLevel(level).fineRatio());

        }
        /////////////////////////////
        // Setup Ascent
        /////////////////////////////
        // Create an instance of Ascent
        ascent::Ascent ascent;
        conduit::Node open_opts;
#ifdef BL_USE_MPI
        open_opts["mpi_comm"] = MPI_Comm_c2f(ParallelDescriptor::Communicator());
#endif

        ascent.open(open_opts);
        ///////////////////////////////////////////////////////////////////
        // Wrap our AMReX Mesh into a Conduit Mesh Blueprint Tree
        ///////////////////////////////////////////////////////////////////
        conduit::Node bp_mesh;
        MultiLevelToBlueprint( amrptr->finestLevel()+1,
                               mfs,
                               var_names,
                               geoms,
                               amrptr->cumTime(),
                               level_steps,
                               ref_ratios,
                               bp_mesh);

        conduit::Node verify_info;
        if(!conduit::blueprint::mesh::verify(bp_mesh,verify_info))
        {
          // verify failed, print error message
          ASCENT_INFO("Error: Mesh Blueprint Verify Failed!");
          // show details of what went awry
          verify_info.print();
        }
        else
        {
          std::cout << " everything A-ok" << std::endl;
        //      verify_info.print();
        }
        // setup actions
        conduit::Node actions;
        ascent.publish(bp_mesh);

        ascent.execute(actions);

    }
    //
    // Write final checkpoint and plotfile.
    //
    if (amrptr->stepOfLastCheckPoint() < amrptr->levelSteps(0))
    {
        amrptr->checkPoint();
    }

    if (amrptr->stepOfLastPlotFile() < amrptr->levelSteps(0))
    {
        amrptr->writePlotFile();
    }

    delete amrptr;

    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_stop = ParallelDescriptor::second() - run_strt;

    ParallelDescriptor::ReduceRealMax(run_stop,IOProc);

    amrex::Print() << "Run time = " << run_stop << std::endl;

    BL_PROFILE_VAR_STOP(pmain);
    BL_PROFILE_REGION_STOP("main()");
    BL_PROFILE_SET_RUN_TIME(run_stop);
    BL_PROFILE_FINALIZE();


    amrex::Finalize();

    return 0;
}
