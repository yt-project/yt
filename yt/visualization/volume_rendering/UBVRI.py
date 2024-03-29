import numpy as np

johnson_filters = {
    "B": {
        "wavelen": np.array(
            [
                3600,
                3650,
                3700,
                3750,
                3800,
                3850,
                3900,
                3950,
                4000,
                4050,
                4100,
                4150,
                4200,
                4250,
                4300,
                4350,
                4400,
                4450,
                4500,
                4550,
                4600,
                4650,
                4700,
                4750,
                4800,
                4850,
                4900,
                4950,
                5000,
                5050,
                5100,
                5150,
                5200,
                5250,
                5300,
                5350,
                5400,
                5450,
                5500,
                5550,
            ],
            dtype="float64",
        ),
        "trans": np.array(
            [
                0.0,
                0.0,
                0.02,
                0.05,
                0.11,
                0.18,
                0.35,
                0.55,
                0.92,
                0.95,
                0.98,
                0.99,
                1.0,
                0.99,
                0.98,
                0.96,
                0.94,
                0.91,
                0.87,
                0.83,
                0.79,
                0.74,
                0.69,
                0.63,
                0.58,
                0.52,
                0.46,
                0.41,
                0.36,
                0.3,
                0.25,
                0.2,
                0.15,
                0.12,
                0.09,
                0.06,
                0.04,
                0.02,
                0.01,
                0.0,
            ],
            dtype="float64",
        ),
    },
    "I": {
        "wavelen": np.array(
            [
                6800,
                6850,
                6900,
                6950,
                7000,
                7050,
                7100,
                7150,
                7200,
                7250,
                7300,
                7350,
                7400,
                7450,
                7500,
                7550,
                7600,
                7650,
                7700,
                7750,
                7800,
                7850,
                7900,
                7950,
                8000,
                8050,
                8100,
                8150,
                8200,
                8250,
                8300,
                8350,
                8400,
                8450,
                8500,
                8550,
                8600,
                8650,
                8700,
                8750,
                8800,
                8850,
                8900,
                8950,
                9000,
                9050,
                9100,
                9150,
                9200,
                9250,
                9300,
                9350,
                9400,
                9450,
                9500,
                9550,
                9600,
                9650,
                9700,
                9750,
                9800,
                9850,
                9900,
                9950,
                10000,
                10050,
                10100,
                10150,
                10200,
                10250,
                10300,
                10350,
                10400,
                10450,
                10500,
                10550,
                10600,
                10650,
                10700,
                10750,
                10800,
                10850,
                10900,
                10950,
                11000,
                11050,
                11100,
                11150,
                11200,
                11250,
                11300,
                11350,
                11400,
                11450,
                11500,
                11550,
                11600,
                11650,
                11700,
                11750,
                11800,
                11850,
            ],
            dtype="float64",
        ),
        "trans": np.array(
            [
                0.0,
                0.0,
                0.01,
                0.01,
                0.01,
                0.04,
                0.08,
                0.13,
                0.17,
                0.21,
                0.26,
                0.3,
                0.36,
                0.4,
                0.44,
                0.49,
                0.56,
                0.6,
                0.65,
                0.72,
                0.76,
                0.84,
                0.9,
                0.93,
                0.96,
                0.97,
                0.97,
                0.98,
                0.98,
                0.99,
                0.99,
                0.99,
                0.99,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.99,
                0.98,
                0.98,
                0.97,
                0.96,
                0.94,
                0.93,
                0.9,
                0.88,
                0.86,
                0.84,
                0.8,
                0.76,
                0.74,
                0.71,
                0.68,
                0.65,
                0.61,
                0.58,
                0.56,
                0.52,
                0.5,
                0.47,
                0.44,
                0.42,
                0.39,
                0.36,
                0.34,
                0.32,
                0.3,
                0.28,
                0.26,
                0.24,
                0.22,
                0.2,
                0.19,
                0.17,
                0.16,
                0.15,
                0.13,
                0.12,
                0.11,
                0.1,
                0.09,
                0.09,
                0.08,
                0.08,
                0.07,
                0.06,
                0.05,
                0.05,
                0.04,
                0.04,
                0.03,
                0.03,
                0.02,
                0.02,
                0.02,
                0.02,
                0.02,
                0.01,
                0.01,
                0.01,
                0.0,
            ],
            dtype="float64",
        ),
    },
    "R": {
        "wavelen": np.array(
            [
                5200,
                5250,
                5300,
                5350,
                5400,
                5450,
                5500,
                5550,
                5600,
                5650,
                5700,
                5750,
                5800,
                5850,
                5900,
                5950,
                6000,
                6050,
                6100,
                6150,
                6200,
                6250,
                6300,
                6350,
                6400,
                6450,
                6500,
                6550,
                6600,
                6650,
                6700,
                6750,
                6800,
                6850,
                6900,
                6950,
                7000,
                7050,
                7100,
                7150,
                7200,
                7250,
                7300,
                7350,
                7400,
                7450,
                7500,
                7550,
                7600,
                7650,
                7700,
                7750,
                7800,
                7850,
                7900,
                7950,
                8000,
                8050,
                8100,
                8150,
                8200,
                8250,
                8300,
                8350,
                8400,
                8450,
                8500,
                8550,
                8600,
                8650,
                8700,
                8750,
                8800,
                8850,
                8900,
                8950,
                9000,
                9050,
                9100,
                9150,
                9200,
                9250,
                9300,
                9350,
                9400,
                9450,
                9500,
            ],
            dtype="float64",
        ),
        "trans": np.array(
            [
                0.0,
                0.01,
                0.02,
                0.04,
                0.06,
                0.11,
                0.18,
                0.23,
                0.28,
                0.34,
                0.4,
                0.46,
                0.5,
                0.55,
                0.6,
                0.64,
                0.69,
                0.71,
                0.74,
                0.77,
                0.79,
                0.81,
                0.84,
                0.86,
                0.88,
                0.9,
                0.91,
                0.92,
                0.94,
                0.95,
                0.96,
                0.97,
                0.98,
                0.99,
                0.99,
                1.0,
                1.0,
                0.99,
                0.98,
                0.96,
                0.94,
                0.92,
                0.9,
                0.88,
                0.85,
                0.83,
                0.8,
                0.77,
                0.73,
                0.7,
                0.66,
                0.62,
                0.57,
                0.53,
                0.49,
                0.45,
                0.42,
                0.39,
                0.36,
                0.34,
                0.31,
                0.27,
                0.22,
                0.19,
                0.17,
                0.15,
                0.13,
                0.12,
                0.11,
                0.1,
                0.08,
                0.07,
                0.06,
                0.06,
                0.05,
                0.04,
                0.04,
                0.03,
                0.03,
                0.02,
                0.02,
                0.02,
                0.01,
                0.01,
                0.01,
                0.01,
                0.0,
            ],
            dtype="float64",
        ),
    },
    "U": {
        "wavelen": np.array(
            [
                3000,
                3050,
                3100,
                3150,
                3200,
                3250,
                3300,
                3350,
                3400,
                3450,
                3500,
                3550,
                3600,
                3650,
                3700,
                3750,
                3800,
                3850,
                3900,
                3950,
                4000,
                4050,
                4100,
                4150,
            ],
            dtype="float64",
        ),
        "trans": np.array(
            [
                0.0,
                0.04,
                0.1,
                0.25,
                0.61,
                0.75,
                0.84,
                0.88,
                0.93,
                0.95,
                0.97,
                0.99,
                1.0,
                0.99,
                0.97,
                0.92,
                0.73,
                0.56,
                0.36,
                0.23,
                0.05,
                0.03,
                0.01,
                0.0,
            ],
            dtype="float64",
        ),
    },
    "V": {
        "wavelen": np.array(
            [
                4600,
                4650,
                4700,
                4750,
                4800,
                4850,
                4900,
                4950,
                5000,
                5050,
                5100,
                5150,
                5200,
                5250,
                5300,
                5350,
                5400,
                5450,
                5500,
                5550,
                5600,
                5650,
                5700,
                5750,
                5800,
                5850,
                5900,
                5950,
                6000,
                6050,
                6100,
                6150,
                6200,
                6250,
                6300,
                6350,
                6400,
                6450,
                6500,
                6550,
                6600,
                6650,
                6700,
                6750,
                6800,
                6850,
                6900,
                6950,
                7000,
                7050,
                7100,
                7150,
                7200,
                7250,
                7300,
                7350,
            ],
            dtype="float64",
        ),
        "trans": np.array(
            [
                0.0,
                0.0,
                0.01,
                0.01,
                0.02,
                0.05,
                0.11,
                0.2,
                0.38,
                0.67,
                0.78,
                0.85,
                0.91,
                0.94,
                0.96,
                0.98,
                0.98,
                0.95,
                0.87,
                0.79,
                0.72,
                0.71,
                0.69,
                0.65,
                0.62,
                0.58,
                0.52,
                0.46,
                0.4,
                0.34,
                0.29,
                0.24,
                0.2,
                0.17,
                0.14,
                0.11,
                0.08,
                0.06,
                0.05,
                0.03,
                0.02,
                0.02,
                0.01,
                0.01,
                0.01,
                0.01,
                0.01,
                0.01,
                0.01,
                0.01,
                0.01,
                0.01,
                0.01,
                0.01,
                0.01,
                0.0,
            ],
            dtype="float64",
        ),
    },
}

for vals in johnson_filters.values():
    wavelen = vals["wavelen"]
    trans = vals["trans"]
    vals["Lchar"] = wavelen[np.argmax(trans)]
