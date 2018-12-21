
def CD(nose_type,M):

    M0=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,
       1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,
       2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,
       3.1,3.2,3.3,3.4,3.5,3.6,2.7,2.8,2.9,4,
       4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5,
       5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6,
       6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7,
       7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.8,7.9,8]


    if nose_type=='3/4_power':
        CD0=[0.02, 0.0109, 0.0087, 0.006, 0.0051, 0.0042, 0.0037, 0.0049, 0.0407, 0.1779,
             0.3182, 0.4498, 0.4870, 0.4566, 0.4170, 0.3858, 0.3672, 0.3517, 0.3379, 0.3253,
             0.3148, 0.3066, 0.2985, 0.2926, 0.2868, 0.2810, 0.2746, 0.2683, 0.2623, 0.2564,
             0.2506, 0.2469, 0.2446, 0.2423, 0.2399, 0.2376, 0.2352, 0.2332, 0.2314, 0.2297,
             0.2282, 0.2268, 0.2255, 0.2251, 0.2247, 0.2241, 0.2228, 0.2216, 0.2204, 0.2196,
             0.2187, 0.2178, 0.2169, 0.2161, 0.2150, 0.2132, 0.2115, 0.2099, 0.2086, 0.2073,
             0.2057, 0.2038, 0.2018, 0.1998, 0.1977, 0.1957, 0.1923, 0.1884, 0.1845, 0.1829,
             0.1817, 0.1805, 0.1805, 0.1809, 0.1814, 0.1822, 0.1830, 0.1837, 0.1852, 0.1868, 0.19]


    if nose_type=='1/2_power':
        CD0=[0.02, 0.0107, 0.0087, 0.0087, 0.0074, 0.0058, 0.006, 0.0085, 0.0254, 0.1247,
             0.2582, 0.3881, 0.4497, 0.4259, 0.391, 0.3587, 0.3466, 0.3362, 0.3269, 0.3182,
             0.31, 0.3035, 0.2976, 0.2928, 0.288, 0.2837, 0.2793, 0.2751, 0.2708, 0.2673,
             0.2641, 0.2608, 0.2585, 0.2568, 0.2552, 0.2533, 0.2512, 0.2492, 0.2476, 0.2465,
             0.2454, 0.2442, 0.2431, 0.242, 0.2409, 0.2403, 0.2398, 0.2393, 0.2381, 0.2367,
             0.2357, 0.2347, 0.2337, 0.2324, 0.2305, 0.2286, 0.2275, 0.2269, 0.2261, 0.225,
             0.2233, 0.2202, 0.2171, 0.2139, 0.2111, 0.2092, 0.2072, 0.2045, 0.2018, 0.199,
             0.1964, 0.1964, 0.1964, 0.1964, 0.1963, 0.1977, 0.1993, 0.2009, 0.2024, 0.2039, 0.2062]

    if nose_type=='Ogive':
        CD0=[0.0213, 0.0149, 0.0115, 0.0087, 0.0063, 0.0049, 0.0049, 0.0082, 0.0178, 0.1028,
             0.2935, 0.4389, 0.535, 0.5398, 0.5238, 0.5003, 0.4792, 0.4581, 0.4357, 0.4131,
             0.3937, 0.3881, 0.3825, 0.3769, 0.3713, 0.3662, 0.3626, 0.359, 0.3554, 0.3517,
             0.3478, 0.3439, 0.34, 0.3349, 0.3296, 0.3246, 0.328, 0.3314, 0.3348, 0.3389,
             0.3433, 0.3469, 0.3473, 0.3476, 0.348, 0.348, 0.348, 0.348, 0.348, 0.348,
             0.3471, 0.3459, 0.3448, 0.3437, 0.3427, 0.3418, 0.3408, 0.3396, 0.3381, 0.3366,
             0.3351, 0.3336, 0.3314, 0.3286, 0.3258, 0.323, 0.3201, 0.3172, 0.3144, 0.3115,
             0.3089, 0.3084, 0.308, 0.3076, 0.3086, 0.3095, 0.3109, 0.3126, 0.3145, 0.3166, 0.3187]


    if nose_type=='HAACK':
        CD0=[0.0213, 0.0121, 0.0098, 0.0083, 0.007, 0.0057, 0.0049, 0.0049, 0.0179, 0.1298,
             0.2891, 0.4969, 0.5447, 0.5255, 0.4983, 0.461, 0.4461, 0.4329, 0.4199, 0.4069,
             0.3951, 0.389, 0.3845, 0.3801, 0.3756, 0.3712, 0.3667, 0.3618, 0.3565, 0.3512,
             0.3495, 0.3488, 0.3465, 0.3442, 0.3419, 0.3396, 0.3373, 0.3361, 0.3351, 0.3341,
             0.3331, 0.3319, 0.3307, 0.3295, 0.3283, 0.327, 0.3258, 0.3245, 0.3243, 0.3242,
             0.3238, 0.3231, 0.3223, 0.3215, 0.3207, 0.3193, 0.3175, 0.3161, 0.3151, 0.3141,
             0.3114, 0.3088, 0.307, 0.3052, 0.303, 0.2997, 0.2963, 0.2933, 0.2911, 0.2889,
             0.2845, 0.2839, 0.2839, 0.2846, 0.2863, 0.2878, 0.2888, 0.2898, 0.2914, 0.2935, 0.2962]


    if nose_type=='Karman':
        CD0=[0.02, 0.0124, 0.0091, 0.0074, 0.0058, 0.0045, 0.0042, 0.0057, 0.0175, 0.1038,
             0.2994, 0.4155, 0.501, 0.4784, 0.443, 0.4059, 0.3897, 0.375, 0.361, 0.3495,
             0.335, 0.3298, 0.3251, 0.3221, 0.3192, 0.3162, 0.3108, 0.3052, 0.298, 0.2944,
             0.2908, 0.2863, 0.2829, 0.2815, 0.28, 0.2781, 0.2763, 0.2738, 0.2711, 0.2691,
             0.2685, 0.2675, 0.2662, 0.2649, 0.2636, 0.2624, 0.2618, 0.2618, 0.2618, 0.2617,
             0.2603, 0.2582, 0.259, 0.2599, 0.2596, 0.257, 0.2548, 0.253, 0.2508, 0.2476,
             0.2449, 0.2426, 0.2411, 0.2392, 0.2364, 0.2335, 0.2304, 0.2277, 0.2251, 0.2231,
             0.2216, 0.22, 0.2189, 0.2189, 0.2207, 0.2223, 0.2226, 0.2233, 0.226, 0.2287, 0.2312]


    if nose_type=='Biconical':
        CD0=[0.02, 0.0118, 0.0082, 0.0052, 0.0024, 0.0025, 0.0037, 0.0043, 0.0214, 0.1188,
              0.2971, 0.4374, 0.4944, 0.4848, 0.4631, 0.4426, 0.4174, 0.3951, 0.3722, 0.3466,
              0.3298, 0.3214, 0.317, 0.3125, 0.3075, 0.3014, 0.2953, 0.2894, 0.2866, 0.2841,
              0.2818, 0.2794, 0.2779, 0.2764, 0.2749, 0.2734, 0.2721, 0.2709, 0.2697, 0.2684,
              0.2679, 0.2675, 0.2672, 0.2667, 0.2658, 0.2648, 0.264, 0.2632, 0.2625, 0.2618,
              0.2611, 0.2604, 0.2598, 0.2592, 0.2586, 0.2581, 0.2579, 0.2579, 0.2578, 0.2574,
              0.2543, 0.2511, 0.248, 0.2449, 0.2417, 0.2386, 0.2355, 0.2324, 0.2296, 0.2276,
              0.2255, 0.224, 0.2244, 0.2248, 0.2255, 0.2264, 0.2274, 0.2284, 0.2306, 0.2334, 0.2362]


    if nose_type=='Conical':
        CD0=[0.02, 0.0143, 0.012, 0.0105, 0.0091, 0.0091, 0.0098, 0.023, 0.1184, 0.2693,
              0.4086, 0.5169, 0.5486, 0.5214, 0.4816, 0.4484, 0.4311, 0.4137, 0.3964, 0.379,
              0.363, 0.3535, 0.3439, 0.3344, 0.3248, 0.3153, 0.3095, 0.3046, 0.2997, 0.2948,
              0.2899, 0.2854, 0.2831, 0.2808, 0.2785, 0.2762, 0.2739, 0.2716, 0.2693, 0.267,
              0.2647, 0.2639, 0.2633, 0.2627, 0.2621, 0.2615, 0.2609, 0.2602, 0.2591, 0.258,
              0.2569, 0.2557, 0.2546, 0.2535, 0.2524, 0.251, 0.2492, 0.2474, 0.2456, 0.2438,
              0.242, 0.2398, 0.2375, 0.2353, 0.233, 0.2307, 0.2285, 0.2261, 0.2234, 0.2206,
              0.2179, 0.2179, 0.2182, 0.2184, 0.2187, 0.2194, 0.2206, 0.2217, 0.2228, 0.2239, 0.225]


    for i in range (1,len(M0)):
        if M0>M[i-1] and M0<M[i]:
            CD=CD0[i-1]+(CD0[i]-CD0[i-1])/(M0[i]-M0[i-1])*(M-M0[i-1])

    return CD
