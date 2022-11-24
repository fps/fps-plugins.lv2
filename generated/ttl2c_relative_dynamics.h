#ifndef plugin_cb_hh
#define plugin_cb_hh

#include <lv2.h>
#include <stdint.h>
 
typedef struct plugin_state plugin_state_t;

typedef struct plugin {
    struct plugin_state *state;
    float *ports[8];
} plugin_t;

enum plugin_port_indices {
    in = 0,
    out = 1,
    t1 = 2,
    t2 = 3,
    strength = 4,
    delay = 5,
    maxratio = 6,
    minratio = 7,
};

typedef struct plugin_port_in {
    float *data;
} plugin_port_in_t;

typedef struct plugin_port_out {
    float *data;
} plugin_port_out_t;

typedef struct plugin_port_t1 {
    float *data;
} plugin_port_t1_t;

typedef struct plugin_port_t2 {
    float *data;
} plugin_port_t2_t;

typedef struct plugin_port_strength {
    float *data;
} plugin_port_strength_t;

typedef struct plugin_port_delay {
    float *data;
} plugin_port_delay_t;

typedef struct plugin_port_maxratio {
    float *data;
} plugin_port_maxratio_t;

typedef struct plugin_port_minratio {
    float *data;
} plugin_port_minratio_t;

     

typedef struct plugin_callbacks
{
    struct plugin* (*const instantiate)(plugin_t *instance, double sample_rate, const char *bundle_path, const LV2_Feature *const *features);
    void (*const connect_port)(plugin_t *instance, uint32_t port, void *data_location);
    void(* const activate)(plugin_t *instance);
    void(* const run)(plugin_t *instance, uint32_t sample_count, const plugin_port_in_t in, const plugin_port_out_t out, const plugin_port_t1_t t1, const plugin_port_t2_t t2, const plugin_port_strength_t strength, const plugin_port_delay_t delay, const plugin_port_maxratio_t maxratio, const plugin_port_minratio_t minratio);
    void(*const deactivate)(plugin_t *instance);
    void(*const cleanup)(plugin_t *instance);
    const void *(*const extension_data)(const char *uri);
} plugin_callbacks_t;

#endif    