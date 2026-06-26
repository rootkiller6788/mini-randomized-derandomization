/**
 * example1.c - USTCONN: Undirected s-t Connectivity in Logspace
 * ==============================================================
 * Demonstrates Reingold's Theorem (USTCONN in L) by testing
 * connectivity on various undirected graphs.
 *
 * Knowledge: L4 (Reingold Theorem), L6 (USTCONN canonical problem).
 * Refs: Reingold (2008) J.ACM 55(4).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "space_derand.h"
#include "logspace_derand.h"

int main(void) {
    printf("=== USTCONN: Logspace Undirected Connectivity ===\n\n");

    /* Example 1: Connected 5-cycle C5 */
    printf("--- C5 Cycle (connected) ---\n");
    size_t edges1[] = {0,1, 1,2, 2,3, 3,4, 4,0};
    bool connected = false;
    if (sd_reingold_ustconn(5, edges1, 5, 0, 3, &connected)) {
        printf("  s=0, t=3: %s\n", connected ? "CONNECTED" : "NOT CONNECTED");
    }

    /* Example 2: Disconnected graph (two components) */
    printf("\n--- Two Components (disconnected) ---\n");
    size_t edges2[] = {0,1, 2,3};
    sd_reingold_ustconn(4, edges2, 2, 0, 3, &connected);
    printf("  s=0, t=3: %s\n", connected ? "CONNECTED" : "NOT CONNECTED");

    /* Example 3: Single edge */
    printf("\n--- Single Edge ---\n");
    size_t edges3[] = {0,1};
    sd_reingold_ustconn(2, edges3, 1, 0, 1, &connected);
    printf("  s=0, t=1: %s\n", connected ? "CONNECTED" : "NOT CONNECTED");

    /* Example 4: Complete graph K5 */
    printf("\n--- K5 Complete Graph ---\n");
    RegularGraph *K5 = sd_complete_graph(5);
    if (K5) {
        printf("  K5: n=%zu, degree=%zu\n", K5->n, K5->d);
        bool is_conn = sd_is_connected(K5);
        size_t diam = sd_shortest_path(K5, 0, 4);
        printf("  Connected: %s\n", is_conn ? "yes" : "no");
        printf("  Diameter: %zu\n", diam);
        sd_regular_graph_free(K5);
    }

    /* Example 5: Component counting */
    printf("\n--- Component Counting ---\n");
    size_t edges5[] = {0,1, 0,2, 3,4, 4,5, 6,7};
    size_t nc = sd_count_components_undirected(8, edges5, 5);
    printf("  Graph with 8 nodes: %zu components\n", nc);

    /* Example 6: Verify undirected property */
    printf("\n--- Verify Undirected ---\n");
    size_t edges6[] = {0,1, 1,0};
    bool is_undirected = sd_verify_undirected(2, edges6, 2);
    printf("  Edge list is undirected: %s\n", is_undirected ? "yes" : "no");

    printf("\n=== Example1 Complete ===\n");
    return 0;
}
