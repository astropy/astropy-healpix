void interpolate_weights(double lon, double lat, int64_t *ring_indices,
                         double *weights, int Nside);

void ring_neighbours(double lon, double lat, int Nside, int64_t *ring1,
                     int64_t *ring2, double *lon1, double *lat1, double *lon2,
                     double *lat2);
