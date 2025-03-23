#include <stdio.h>
#include "calculator.h"
#include "advanced_calculator.h"

int main() {
    printf("Calculator Library Test\n");
    printf("----------------------\n\n");
    
    // 测试基本算术操作
    printf("Basic Arithmetic:\n");
    printf("10 + 5 = %.2f\n", add(10, 5));
    printf("10 - 5 = %.2f\n", subtract(10, 5));
    printf("10 * 5 = %.2f\n", multiply(10, 5));
    printf("10 / 5 = %.2f\n", divide(10, 5));
    printf("\n");
    
    // 测试高级数学函数
    printf("Advanced Math:\n");
    printf("2^10 = %.2f\n", power(2, 10));
    printf("sqrt(16) = %.2f\n", square_root(16));
    printf("factorial(5) = %d\n", factorial(5));
    printf("abs(-7.5) = %.2f\n", absolute(-7.5));
    printf("\n");
    
    // 测试三角函数
    printf("Trigonometry:\n");
    printf("sin(30°) = %.4f\n", sine(30));
    printf("cos(60°) = %.4f\n", cosine(60));
    printf("tan(45°) = %.4f\n", tangent(45));
    printf("\n");
    
    // 测试对数函数
    printf("Logarithmic Functions:\n");
    printf("ln(10) = %.4f\n", natural_log(10));
    printf("log10(100) = %.4f\n", log_base10(100));
    printf("\n");
    
    // 测试统计函数
    printf("Statistics:\n");
    double values[] = {1, 2, 3, 4, 5};
    int size = sizeof(values) / sizeof(values[0]);
    printf("Mean of [1,2,3,4,5] = %.2f\n", calculate_mean(values, size));
    printf("Variance of [1,2,3,4,5] = %.2f\n", calculate_variance(values, size));
    printf("Standard Deviation of [1,2,3,4,5] = %.2f\n", calculate_std_deviation(values, size));
    printf("\n");
    
    // 测试转换函数
    printf("Conversions:\n");
    printf("25°C in Fahrenheit = %.2f°F\n", celsius_to_fahrenheit(25));
    printf("77°F in Celsius = %.2f°C\n", fahrenheit_to_celsius(77));
    printf("180° in radians = %.4f rad\n", degrees_to_radians(180));
    printf("3.14159 rad in degrees = %.2f°\n", radians_to_degrees(3.14159));
    printf("\n");
    
    // 测试金融函数
    printf("Financial Calculations:\n");
    printf("Simple interest (1000, 5%%, 2 years) = %.2f\n", calculate_simple_interest(1000, 0.05, 2));
    printf("Compound interest (1000, 5%%, 2 years, 12 compoundings) = %.2f\n", 
           calculate_compound_interest(1000, 0.05, 2, 12));
    printf("Monthly payment on 10000 loan at 6%% for 36 months = %.2f\n", 
           calculate_loan_payment(10000, 0.06, 36));
    printf("\n");
    
    // 测试位操作函数
    printf("Bit Operations:\n");
    printf("Number of bits in 255 = %u\n", bit_count(255));
    printf("Bit reversal of 0x0F00 = 0x%X\n", bit_reverse(0x0F00));
    printf("Is 64 a power of 2? %s\n", is_power_of_two(64) ? "Yes" : "No");
    printf("Is 63 a power of 2? %s\n", is_power_of_two(63) ? "Yes" : "No");
    printf("\n");
    
    // 测试矩阵操作
    printf("Matrix Operations:\n");
    Matrix* m1 = matrix_create(2, 2);
    if (m1) {
        m1->data[0][0] = 1; m1->data[0][1] = 2;
        m1->data[1][0] = 3; m1->data[1][1] = 4;
        
        printf("Determinant of [[1,2],[3,4]] = %.2f\n", matrix_determinant(m1));
        
        Matrix* transpose = matrix_transpose(m1);
        if (transpose) {
            printf("Transpose of [[1,2],[3,4]] = [[%.0f,%.0f],[%.0f,%.0f]]\n",
                   transpose->data[0][0], transpose->data[0][1],
                   transpose->data[1][0], transpose->data[1][1]);
            matrix_free(transpose);
        }
        
        matrix_free(m1);
    }
    printf("\n");
    
    // 测试复数操作
    printf("Complex Number Operations:\n");
    Complex a = {3, 4}; // 3 + 4i
    Complex b = {1, 2}; // 1 + 2i
    
    Complex sum = complex_add(a, b);
    Complex product = complex_multiply(a, b);
    
    printf("(%g + %gi) + (%g + %gi) = %g + %gi\n", a.real, a.imag, b.real, b.imag, sum.real, sum.imag);
    printf("(%g + %gi) * (%g + %gi) = %g + %gi\n", a.real, a.imag, b.real, b.imag, product.real, product.imag);
    printf("|%g + %gi| = %g\n", a.real, a.imag, complex_modulus(a));
    printf("\n");
    
    return 0;
} 