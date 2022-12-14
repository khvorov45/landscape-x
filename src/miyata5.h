int  locpenaltym = -1440;
int  exgpm = +0; /* != 0 nisuruto kowareru. exgp ha constants.c de kurikomu */
char locaminom[] = "ARNDCQEGHILKMFPSTWYVBZX.-J";
char locgrpm[] =
    {
        0,
        3,
        2,
        2,
        5,
        2,
        2,
        0,
        3,
        1,
        1,
        3,
        1,
        4,
        0,
        0,
        0,
        4,
        4,
        1,
        2,
        2,
        6,
        6,
        6,
        1,
};
int locn_dism[26][26] =
    {
        {
            600,
            -235,
            91,
            -78,
            202,
            51,
            -103,
            340,
            -21,
            -169,
            -189,
            -246,
            -92,
            -323,
            582,
            454,
            342,
            -400,
            -309,
            71,
            7,
            -26,
            -15,
            -400,
            0,
            -1400,
        },

        {
            -235,
            600,
            17,
            -69,
            -275,
            277,
            185,
            -400,
            365,
            -112,
            -149,
            485,
            -55,
            -106,
            -229,
            -183,
            20,
            -178,
            22,
            -95,
            -26,
            231,
            -15,
            -400,
            0,
            -1400,
        },

        {
            91,
            17,
            600,
            414,
            -209,
            317,
            357,
            39,
            231,
            -363,
            -398,
            74,
            -280,
            -400,
            85,
            225,
            200,
            -400,
            -378,
            -189,
            507,
            337,
            -15,
            -400,
            0,
            -1400,
        },

        {
            -78,
            -69,
            414,
            600,
            -395,
            179,
            342,
            -78,
            108,
            -400,
            -400,
            14,
            -400,
            -400,
            -86,
            65,
            14,
            -400,
            -400,
            -372,
            507,
            261,
            -15,
            -400,
            0,
            -1400,
        },

        {
            202,
            -275,
            -209,
            -395,
            600,
            -109,
            -332,
            -35,
            -132,
            134,
            128,
            -335,
            182,
            -40,
            220,
            74,
            185,
            -355,
            -81,
            354,
            -302,
            -220,
            -15,
            -400,
            0,
            -1400,
        },

        {
            51,
            277,
            317,
            179,
            -109,
            600,
            360,
            -109,
            508,
            -135,
            -172,
            297,
            -58,
            -203,
            51,
            128,
            280,
            -378,
            -109,
            -9,
            248,
            480,
            -15,
            -400,
            0,
            -1400,
        },

        {
            -103,
            185,
            357,
            342,
            -332,
            360,
            600,
            -195,
            325,
            -369,
            -400,
            274,
            -295,
            -400,
            -109,
            11,
            77,
            -400,
            -321,
            -249,
            350,
            480,
            -15,
            -400,
            0,
            -1400,
        },

        {
            340,
            -400,
            39,
            -78,
            -35,
            -109,
            -195,
            600,
            -195,
            -400,
            -400,
            -400,
            -355,
            -400,
            322,
            357,
            114,
            -400,
            -400,
            -189,
            -19,
            -152,
            -15,
            -400,
            0,
            -1400,
        },

        {
            -21,
            365,
            231,
            108,
            -132,
            508,
            325,
            -195,
            600,
            -100,
            -141,
            374,
            -26,
            -152,
            -15,
            45,
            222,
            -303,
            -49,
            -3,
            169,
            417,
            -15,
            -400,
            0,
            -1400,
        },

        {
            -169,
            -112,
            -363,
            -400,
            134,
            -135,
            -369,
            -400,
            -100,
            600,
            560,
            -212,
            517,
            425,
            -149,
            -243,
            -12,
            108,
            354,
            357,
            -400,
            -252,
            -15,
            -400,
            0,
            -1400,
        },

        {
            -189,
            -149,
            -398,
            -400,
            128,
            -172,
            -400,
            -400,
            -141,
            560,
            600,
            -252,
            482,
            420,
            -172,
            -269,
            -43,
            105,
            331,
            340,
            -400,
            -290,
            -15,
            -400,
            0,
            -1400,
        },

        {
            -246,
            485,
            74,
            14,
            -335,
            297,
            274,
            -400,
            374,
            -212,
            -252,
            600,
            -152,
            -215,
            -240,
            -175,
            -1,
            -289,
            -92,
            -172,
            44,
            285,
            -15,
            -400,
            0,
            -1400,
        },

        {
            -92,
            -55,
            -280,
            -400,
            182,
            -58,
            -295,
            -355,
            -26,
            517,
            482,
            -152,
            600,
            365,
            -75,
            -163,
            68,
            59,
            334,
            422,
            -368,
            -176,
            -15,
            -400,
            0,
            -1400,
        },

        {
            -323,
            -106,
            -400,
            -400,
            -40,
            -203,
            -400,
            -400,
            -152,
            425,
            420,
            -215,
            365,
            600,
            -306,
            -386,
            -143,
            282,
            462,
            191,
            -400,
            -315,
            -15,
            -400,
            0,
            -1400,
        },

        {
            582,
            -229,
            85,
            -86,
            220,
            51,
            -109,
            322,
            -15,
            -149,
            -172,
            -240,
            -75,
            -306,
            600,
            440,
            351,
            -400,
            -292,
            88,
            0,
            -29,
            -15,
            -400,
            0,
            -1400,
        },

        {
            454,
            -183,
            225,
            65,
            74,
            128,
            11,
            357,
            45,
            -243,
            -269,
            -175,
            -163,
            -386,
            440,
            600,
            345,
            -400,
            -352,
            -15,
            145,
            70,
            -15,
            -400,
            0,
            -1400,
        },

        {
            342,
            20,
            200,
            14,
            185,
            280,
            77,
            114,
            222,
            -12,
            -43,
            -1,
            68,
            -143,
            351,
            345,
            600,
            -400,
            -100,
            194,
            107,
            178,
            -15,
            -400,
            0,
            -1400,
        },

        {
            -400,
            -178,
            -400,
            -400,
            -355,
            -378,
            -400,
            -400,
            -303,
            108,
            105,
            -289,
            59,
            282,
            -400,
            -400,
            -400,
            600,
            297,
            -118,
            -400,
            -400,
            -15,
            -400,
            0,
            -1400,
        },

        {
            -309,
            22,
            -378,
            -400,
            -81,
            -109,
            -321,
            -400,
            -49,
            354,
            331,
            -92,
            334,
            462,
            -292,
            -352,
            -100,
            297,
            600,
            165,
            -400,
            -215,
            -15,
            -400,
            0,
            -1400,
        },

        {
            71,
            -95,
            -189,
            -372,
            354,
            -9,
            -249,
            -189,
            -3,
            357,
            340,
            -172,
            422,
            191,
            88,
            -15,
            194,
            -118,
            165,
            600,
            -280,
            -129,
            -15,
            -400,
            0,
            -1400,
        },

        {
            7,
            -26,
            507,
            507,
            -302,
            248,
            350,
            -19,
            169,
            -400,
            -400,
            44,
            -368,
            -400,
            0,
            145,
            107,
            -400,
            -400,
            -280,
            507,
            299,
            -400,
            -400,
            0,
            -1400,
        },

        {
            -26,
            231,
            337,
            261,
            -220,
            480,
            480,
            -152,
            417,
            -252,
            -290,
            285,
            -176,
            -315,
            -29,
            70,
            178,
            -400,
            -215,
            -129,
            299,
            480,
            -400,
            -400,
            0,
            -1400,
        },

        {
            -15,
            -15,
            -15,
            -15,
            -15,
            -15,
            -15,
            -15,
            -15,
            -15,
            -15,
            -15,
            -15,
            -15,
            -15,
            -15,
            -15,
            -15,
            -15,
            -15,
            -400,
            -400,
            -400,
            -400,
            0,
            -1400,
        },

        {
            -400,
            -400,
            -400,
            -400,
            -400,
            -400,
            -400,
            -400,
            -400,
            -400,
            -400,
            -400,
            -400,
            -400,
            -400,
            -400,
            -400,
            -400,
            -400,
            -400,
            -400,
            -400,
            -400,
            -400,
            0,
            -1400,
        },

        {
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        },

        {
            -1400,
            -1400,
            -1400,
            -1400,
            -1400,
            -1400,
            -1400,
            -1400,
            -1400,
            -1400,
            -1400,
            -1400,
            -1400,
            -1400,
            -1400,
            -1400,
            -1400,
            -1400,
            -1400,
            -1400,
            -1400,
            -1400,
            -1400,
            -1400,
            0,
            1600,
        },
};
