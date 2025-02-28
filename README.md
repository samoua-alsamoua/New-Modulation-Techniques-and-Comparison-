# **New Modulation Techniques and Comparison**

## **Overview**
This project implements and compares various modulation techniques used in modern communication systems. The goal is to evaluate the performance of traditional modulation schemes (e.g., BPSK, QPSK) and advanced techniques (e.g., 16-QAM, 64-QAM, OFDM) under different Signal-to-Noise Ratio (SNR) conditions. The performance metrics include **Bit Error Rate (BER)** and **spectral efficiency**.

The project also demonstrates the implementation of **Orthogonal Frequency Division Multiplexing (OFDM)**, a widely used technique in wireless communication systems like Wi-Fi, LTE, and 5G.

---

## **Features**
- **Modulation Techniques**: Implements BPSK, QPSK, 16-QAM, 64-QAM, and OFDM.
- **Performance Metrics**: Compares BER and spectral efficiency across modulation schemes.
- **Visualization**: Includes plots for BER vs SNR and spectral efficiency comparison.
- **Scalability**: Easily extendable to test other modulation schemes or channel conditions (e.g., multipath fading).

---

## **Requirements**
- MATLAB toolboxes: **Communications Toolbox** and **Signal Processing Toolbox**.

---

## **Code Structure**
The code is organized into the following sections:

### **1. Parameters**
Defines simulation parameters such as:
- Modulation orders (`M_qam` for QAM).
- SNR range (`snr_range`).
- Number of bits to simulate (`num_bits`).

### **2. Modulation Implementation**
Implements the following modulation techniques:
- **BPSK**: Binary Phase Shift Keying.
- **QPSK**: Quadrature Phase Shift Keying.
- **16-QAM**: 16-level Quadrature Amplitude Modulation.
- **64-QAM**: 64-level Quadrature Amplitude Modulation.
- **OFDM**: Orthogonal Frequency Division Multiplexing with QPSK modulation on subcarriers.

### **3. Performance Metrics**
Calculates:
- **BER**: Bit Error Rate under varying SNR conditions using MATLAB's `biterr` function.
- **Spectral Efficiency**: Bits per symbol for each modulation scheme.

### **4. Visualization**
Generates:
- BER vs SNR plot.
- Spectral efficiency bar chart.

---

## **View Results**:
   - The script will generate two plots:
     1. **BER vs SNR**: A semilog plot comparing the BER performance of different modulation techniques.
     2. **Spectral Efficiency**: A bar chart comparing the spectral efficiency of each modulation scheme.
   - Additionally, computational complexity details will be printed in the MATLAB Command Window.

### **BER vs SNR Plot**
The plot shows the BER performance of different modulation techniques under varying SNR conditions. Lower-order modulations (e.g., BPSK, QPSK) perform better in noisy environments, while higher-order modulations (e.g., 64-QAM) require higher SNR for reliable communication.

### **Spectral Efficiency Comparison**
The bar chart compares the spectral efficiency of each modulation scheme. Higher-order modulations (e.g., 64-QAM) achieve higher spectral efficiency but are less robust to noise.

---

## **Acknowledgments**
- MATLAB Documentation: [Communications Toolbox](https://www.mathworks.com/help/comm/)
- Inspiration from modern communication systems like Wi-Fi, LTE, and 5G.

---

## Contact
- **Author**: Samoua Alsamoua
- **GitHub**: [samoua-alsamoua](https://github.com/samoua-alsamoua)
- **Website**: [saalsamoua.github.io](https://samoua-alsamoua.github.io/saalsamoua/)

```

---
