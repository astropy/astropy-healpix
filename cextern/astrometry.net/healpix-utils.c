/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include "bl.h"
#include "healpix.h"
#include "mathutil.h"
#include "starutil.h"
#include <stdio.h>

ll* healpix_region_search(int seed, ll* seeds, int Nside,
                          ll* accepted, ll* rejected,
                          int (*accept)(int hp, void* token),
                          void* token, int depth) {
    ll* frontier;
    anbool allocd_rej = FALSE;
    int d;

    if (!accepted)
        accepted = ll_new(256);
    if (!rejected) {
        rejected = ll_new(256);
        allocd_rej = TRUE;
    }

    if (seeds)
        //frontier = seeds;
        frontier = ll_dupe(seeds);
    else {
        frontier = ll_new(256);
        ll_append(frontier, seed);
    }

    for (d=0; !depth || d<depth; d++) {
        int j, N;
        N = ll_size(frontier);
        if (N == 0)
            break;
        for (j=0; j<N; j++) {
            int64_t hp;
            int i;
            int64_t neigh[8];
            hp = ll_get(frontier, j);
            healpixl_get_neighbours(hp, neigh, Nside);
            for (i=0; i<8; i++) {
                if (neigh[i] < 0)
                    continue;
                if (ll_contains(frontier, neigh[i]))
                    continue;
                if (ll_contains(rejected, neigh[i]))
                    continue;
                if (ll_contains(accepted, neigh[i]))
                    continue;
                if (accept(neigh[i], token)) {
                    ll_append(accepted, neigh[i]);
                    ll_append(frontier, neigh[i]);
                } else
                    ll_append(rejected, neigh[i]);
            }
        }
        ll_remove_index_range(frontier, 0, N);
    }

    ll_free(frontier);
    if (allocd_rej)
        ll_free(rejected);
    return accepted;
}


static ll* hp_rangesearch(const double* xyz, double radius, int Nside, ll* hps, anbool approx) {
    int64_t hp;
    double hprad = arcmin2dist(healpix_side_length_arcmin(Nside)) * sqrt(2);
    ll* frontier = ll_new(256);
    ll* bad = ll_new(256);
    if (!hps)
        hps = ll_new(256);

    hp = xyzarrtohealpixl(xyz, Nside);
    ll_append(frontier, hp);
    ll_append(hps, hp);
    while (ll_size(frontier)) {
        int64_t neighbours[8];
        int64_t i;
        hp = ll_pop(frontier);
        healpixl_get_neighbours(hp, neighbours, Nside);
        for (i=0; i<8; i++) {
            anbool tst;
            double nxyz[3];
            if (neighbours[i] < 0)
                continue;
            if (ll_contains(frontier, neighbours[i]))
                continue;
            if (ll_contains(bad, neighbours[i]))
                continue;
            if (ll_contains(hps, neighbours[i]))
                continue;
            if (approx) {
                healpixl_to_xyzarr(neighbours[i], Nside, 0.5, 0.5, nxyz);
                tst = (sqrt(distsq(xyz, nxyz, 3)) - hprad <= radius);
            } else {
                tst = healpixl_within_range_of_xyz(neighbours[i], Nside, xyz, radius);
            }
            if (tst) {
                // in range!
                ll_append(frontier, neighbours[i]);
                ll_append(hps, neighbours[i]);
            } else
                ll_append(bad, neighbours[i]);
        }
    }

    ll_free(bad);
    ll_free(frontier);

    return hps;
}

ll* healpix_rangesearch_xyz_approx(const double* xyz, double radius, int Nside, ll* hps) {
    return hp_rangesearch(xyz, radius, Nside, hps, TRUE);
}

ll* healpix_rangesearch_xyz(const double* xyz, double radius, int Nside, ll* hps) {
    return hp_rangesearch(xyz, radius, Nside, hps, FALSE);
}

ll* healpix_rangesearch_radec_approx(double ra, double dec, double radius, int Nside, ll* hps) {
    double xyz[3];
    radecdeg2xyzarr(ra, dec, xyz);
    return hp_rangesearch(xyz, radius, Nside, hps, TRUE);
}

ll* healpix_rangesearch_radec(double ra, double dec, double radius, int Nside, ll* hps) {
    double xyz[3];
    radecdeg2xyzarr(ra, dec, xyz);
    return hp_rangesearch(xyz, radius, Nside, hps, FALSE);
}

// The following is a version of the healpix_rangesearch_radec function
// that works with standard C types to make interfacing with Python easier.
int64_t healpix_rangesearch_radec_simple(double ra, double dec, double radius,
                                         int Nside, int approx, int64_t **indices) {
    double xyz[3];
    ll* hps = ll_new(256);
    int64_t n_pixels;
    radecdeg2xyzarr(ra, dec, xyz);
    hp_rangesearch(xyz, radius, Nside, hps, (anbool)approx);
    n_pixels = ll_size(hps);
    *indices = (int64_t *)malloc((size_t)ll_size(hps) * sizeof(int64_t));
    if (*indices) {
        ll_copy(hps, 0, hps->N, *indices);
    } else {
        fprintf(stderr, "malloc failed\n");
        n_pixels = -1;
    }
    ll_free(hps);
    return n_pixels;
}
