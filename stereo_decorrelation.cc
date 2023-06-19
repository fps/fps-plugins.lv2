#include <atomic>
#include <algorithm>

#include <lv2.h>
#include <lv2/worker/worker.h>
#include <lv2/urid/urid.h>
#include <lv2/atom/atom.h>

#include <FFTConvolver/FFTConvolver.h>

#include "common.h"

#include "stereo_decorrelation.h"

#define STEREO_DECORRELATION_URI FPS_PLUGINS_BASE_URI "/stereo_decorrelation"

// #define STEREO_DECORRELATION_LOG(x) { std::cerr << "stereo_decorrelation: " << x; }
#define STEREO_DECORRELATION_LOG(x) { }

#define STEREO_DECORRELATION_BLOCK_SIZE 32

struct plugin
{
    std::atomic<bool> m_working;
    LV2_Worker_Schedule m_worker_schedule;
    LV2_URID_Map m_urid_map;
    float m_sample_rate;
    
    std::vector<float *> m_ports;

    stereo_decorrelation m_stereo_decorrelation;

    int m_previous_method;
    bool m_previous_whiten;
    bool m_previous_sum_to_mono;
    float m_previous_decay;
    int m_previous_seed_left;
    int m_previous_seed_right;

    std::vector<float> m_input_buffer_left;
    std::vector<float> m_input_buffer_right;
    std::vector<float> m_output_buffer_left;
    std::vector<float> m_output_buffer_right;
    
    fftconvolver::FFTConvolver m_left_convolver;
    fftconvolver::FFTConvolver m_right_convolver;

    plugin (float sample_rate) :
        m_working (false),
        m_sample_rate (sample_rate),
        m_ports (11, 0),
        m_stereo_decorrelation (0.1 * sample_rate),

        m_previous_method (0),
        m_previous_whiten (true),
        m_previous_sum_to_mono (false),
        m_previous_decay (0),
        m_previous_seed_left (11),
        m_previous_seed_right (12),

        m_input_buffer_left (STEREO_DECORRELATION_BLOCK_SIZE, 0),
        m_input_buffer_right (STEREO_DECORRELATION_BLOCK_SIZE, 0),
        m_output_buffer_left (STEREO_DECORRELATION_BLOCK_SIZE, 0),
        m_output_buffer_right (STEREO_DECORRELATION_BLOCK_SIZE, 0)
    {
        STEREO_DECORRELATION_LOG("plugin ()\n");
        init (m_previous_method, m_previous_whiten, m_previous_sum_to_mono, m_previous_decay, m_previous_seed_left, m_previous_seed_right);
    }
    
    void init
    (
        int method, bool whiten, bool sum_to_mono,
        float decay, int seed_left, int seed_right
    )
    {
        m_working = true;

        m_previous_method = method;
        m_previous_whiten = whiten;
        m_previous_sum_to_mono = sum_to_mono;
        m_previous_decay = decay;
        m_previous_seed_left = seed_left;
        m_previous_seed_right = seed_right;

        STEREO_DECORRELATION_LOG("plugin::init ()\n");
        m_stereo_decorrelation.init (m_sample_rate*decay, seed_left, seed_right, whiten, sum_to_mono);
        
        m_left_convolver.init 
        (
            STEREO_DECORRELATION_BLOCK_SIZE, 
            &m_stereo_decorrelation.m_left_response[0], 
            m_stereo_decorrelation.m_left_response.size ()
        );
        
        m_right_convolver.init 
        (
            STEREO_DECORRELATION_BLOCK_SIZE, 
            &m_stereo_decorrelation.m_right_response[0], 
            m_stereo_decorrelation.m_left_response.size ()
        );
    }
};

LV2_Handle instantiate
(
    const LV2_Descriptor *descriptor,
    double sample_rate,
    const char *bundle_path,
    const LV2_Feature *const *features
)
{
    STEREO_DECORRELATION_LOG("instantiate\n")
    const LV2_Feature * feature = 0;
    bool worker_found = false;
    LV2_Worker_Schedule worker_schedule;
    LV2_URID_Map urid_map;

    bool urid_map_found = false;
    STEREO_DECORRELATION_LOG("instantiate2\n")
    while ((feature = *(++features)) != 0)
    {
        STEREO_DECORRELATION_LOG(feature->URI << "\n")
        if (std::string (feature->URI) == LV2_WORKER__schedule)
        {
            worker_found = true;
            worker_schedule = *((LV2_Worker_Schedule*)feature->data);
        }

        if (std::string (feature->URI) == LV2_URID__map)
        {
            urid_map_found = true;
            urid_map = *((LV2_URID_Map*)feature->data);
        }
    }
    STEREO_DECORRELATION_LOG("instantiate3\n")

    if (false == worker_found) return 0;
    if (false == urid_map_found) return 0;

            
    STEREO_DECORRELATION_LOG("instantiate4\n")

    plugin *instance = new plugin (sample_rate);
    STEREO_DECORRELATION_LOG("instantiate5\n")

    instance->m_worker_schedule = worker_schedule;
    instance->m_urid_map = urid_map;
    STEREO_DECORRELATION_LOG("done. returning\n")
    return (LV2_Handle)instance;
}


static void cleanup(LV2_Handle instance)
{
    STEREO_DECORRELATION_LOG("cleanup\n")
    // fftw_cleanup ();
    delete ((plugin*)instance);
}

static void connect_port (LV2_Handle instance, uint32_t port, void *data_location)
{
        plugin &p = *((plugin*)instance);
        p.m_ports[port] = (float*)data_location;
}

static void run
(
    LV2_Handle instance,
    uint32_t sample_count
)
{
    //STEREO_DECORRELATION_LOG("run")
    plugin &the_plugin = *((plugin*)instance);

    // Audio ports
    const float *inl           =  the_plugin.m_ports[0];
    const float *inr           =  the_plugin.m_ports[1];
    float       *outl          =  the_plugin.m_ports[2];
    float       *outr          =  the_plugin.m_ports[3];

    // Control eports
    const float &amount        = *the_plugin.m_ports[4];
    const float &dry_amount    = *the_plugin.m_ports[5];
    const float &method        = *the_plugin.m_ports[6];
    const float &sum_to_mono   = *the_plugin.m_ports[7];
    const float &decay         = *the_plugin.m_ports[8];
    const float &seed_left     = *the_plugin.m_ports[9];
    const float &seed_right    = *the_plugin.m_ports[10];

    float data[5] = { decay, seed_left, seed_right, whiten, sum_to_mono };

    if (the_plugin.m_working) 
    {
        for (size_t index = 0; index < sample_count; ++index)
        {
            outl[index] = 0;
            outr[index] = 0;
        }
        return;
    }
    
    if
    (
        the_plugin.m_previous_decay != decay ||
        the_plugin.m_previous_seed_left != seed_left ||
        the_plugin.m_previous_seed_right != seed_right ||
        the_plugin.m_pre
    )
    {
        STEREO_DECORRELATION_LOG ("scheduling work\n");
        
        the_plugin.m_working = true;

        the_plugin.m_worker_schedule.schedule_work
        (
            the_plugin.m_worker_schedule.handle,
            5 * sizeof (float),
            &data
        );
        
        for (size_t index = 0; index < sample_count; ++index)
        {
            outl[index] = 0;
            outr[index] = 0;
        }
        return;
    }
    
    // STEREO_DECORRELATION_LOG("not working\n");
    
    size_t processed_samples = 0;
    size_t remaining_samples = sample_count;
    
    while (remaining_samples != 0)
    {
        // STEREO_DECORRELATION_LOG("remaining: " << remaining_samples << "\n")
        
        size_t samples_to_process = (size_t)std::min ((size_t)STEREO_DECORRELATION_BLOCK_SIZE, remaining_samples);
        
        // STEREO_DECORRELATION_LOG("to process: " << samples_to_process << "\n")
        
        std::copy 
        (
            inl + processed_samples, 
            inl + processed_samples + samples_to_process, 
            the_plugin.m_input_buffer_left.begin ()
        );

        std::copy 
        (
            inr + processed_samples, 
            inr + processed_samples + samples_to_process, 
            the_plugin.m_input_buffer_right.begin ()
        );

        the_plugin.m_left_convolver.process 
        (
            &the_plugin.m_input_buffer_left[0], 
            &the_plugin.m_output_buffer_left[0], 
            samples_to_process
        );
        
        the_plugin.m_right_convolver.process 
        (
            &the_plugin.m_input_buffer_right[0], 
            &the_plugin.m_output_buffer_right[0], 
            samples_to_process
        );
        
        const float amount_gain_factor = pow(10, amount/20);
        const float dry_amount_gain_factor = pow(10, dry_amount/20);
        for (size_t index = 0; index < samples_to_process; ++index)
        {
            outl[index + processed_samples] = amount_gain_factor * the_plugin.m_output_buffer_left[index] + dry_amount_gain_factor * inl[index + processed_samples];
            
            outr[index + processed_samples] = amount_gain_factor * the_plugin.m_output_buffer_right[index] + dry_amount_gain_factor * inr[index + processed_samples];
        }
        
        processed_samples += samples_to_process;
        remaining_samples -= samples_to_process;
    }
}


LV2_Worker_Status work
(
    LV2_Handle instance,
    LV2_Worker_Respond_Function respond,
    LV2_Worker_Respond_Handle handle,
    uint32_t size,
    const void *data
)
{
    STEREO_DECORRELATION_LOG("work: " << size << "\n")

    if (size != 2 * sizeof (float))
    {
        std::cerr << "stereo_decorrelation: Bad data!\n";
        return LV2_WORKER_ERR_UNKNOWN;
    }
    
    plugin &the_plugin = *((plugin*)instance);

    float decay = ((float*)data)[0];
    float seed = ((float*)data)[1];
    
    STEREO_DECORRELATION_LOG("decay: " << decay << " seed: " << seed << "\n");
    
    the_plugin.init (decay, seed);
    
    the_plugin.m_working = false;

    return LV2_WORKER_SUCCESS;
}

LV2_Worker_Status work_response
(
    LV2_Handle instance,
    uint32_t size,
    const void *body
)
{
    return LV2_WORKER_SUCCESS;
}

LV2_Worker_Status end_run(LV2_Handle instance)
{
    return LV2_WORKER_SUCCESS;
}

static LV2_Worker_Interface worker_interface =
{
    .work = work,
    .work_response = work_response,
    .end_run = end_run,
};

static const void *extension_data (const char *uri)
{
    STEREO_DECORRELATION_LOG("get extension_data. URI: " << uri << "\n")

    if (std::string(uri) == LV2_WORKER__interface)
    {
        STEREO_DECORRELATION_LOG("get extension data for worker\n")
        return &worker_interface;
    }

    return 0;
}

static LV2_Descriptor plugin_descriptor = {
    STEREO_DECORRELATION_URI,
    instantiate,
    connect_port,
    0,
    run,
    0,
    cleanup,
    extension_data
};

LV2_SYMBOL_EXPORT const LV2_Descriptor* lv2_descriptor (uint32_t index) {
    if (0 == index) {
          return &plugin_descriptor;
    } else {
        return NULL;
    }
}

