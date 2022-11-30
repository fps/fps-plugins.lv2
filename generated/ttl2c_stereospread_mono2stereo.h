#ifndef plugin_cb_hh
#define plugin_cb_hh

#include <lv2.h>
#include <lv2/log/logger.h>
#include <lv2/core/lv2_util.h>
#include <stdint.h>
 
typedef struct plugin_state plugin_state_t;

typedef struct {
    plugin_state_t *state;
    void *ports[5];
    LV2_URID_Map *map;
    LV2_Log_Logger logger;
} plugin_t;

enum plugin_port_indices {
    in = 0,
    outl = 1,
    outr = 2,
    drywet = 3,
    length = 4,
};

typedef struct {
    float const *data;
} plugin_port_in_t;

typedef struct {
    float *data;
} plugin_port_outl_t;

typedef struct {
    float *data;
} plugin_port_outr_t;

typedef struct {
    float const data;
} plugin_port_drywet_t;

typedef struct {
    float const data;
} plugin_port_length_t;

typedef struct {
    plugin_t* (*const instantiate)(plugin_t *instance, double sample_rate, const char *bundle_path, const LV2_Feature *const *features);
    void (*const connect_port)(plugin_t *instance, uint32_t port, void *data_location);
    void (*const activate)(plugin_t *instance);
    void (*const run)(plugin_t *instance, uint32_t sample_count, const plugin_port_in_t in, const plugin_port_outl_t outl, const plugin_port_outr_t outr, const plugin_port_drywet_t drywet, const plugin_port_length_t length);
    void (*const deactivate)(plugin_t *instance);
    void (*const cleanup)(plugin_t *instance);
    const void *(*const extension_data)(const char *uri);
} plugin_callbacks_t;

#endif    