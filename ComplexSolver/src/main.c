#include <ti/getcsc.h>
#include <ti/screen.h>
#include <ti/real.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* --- GLOBALS (Prevents Stack Overflow) --- */
#define MAX_N 5
typedef struct
{
    float r;
    float i;
} complex_t;

complex_t A[MAX_N][MAX_N + 1]; // Main Matrix
complex_t M[MAX_N][MAX_N + 1]; // Solver Temp Matrix
complex_t X[MAX_N];            // Solution Vector
char input_buf[60];            // Global Input Buffer
char num_buf[32];              // Global Number Parser Buffer
char fmt_buf[60];              // Increased buffer for long strings
int n = 2;
int cur_r = 0;
int cur_c = 0;

/* --- MATH HELPERS --- */
complex_t c_sub(complex_t a, complex_t b) { return (complex_t){a.r - b.r, a.i - b.i}; }
complex_t c_add(complex_t a, complex_t b) { return (complex_t){a.r + b.r, a.i + b.i}; }
complex_t c_mul(complex_t a, complex_t b) { return (complex_t){a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r}; }

complex_t c_div(complex_t a, complex_t b)
{
    float d = b.r * b.r + b.i * b.i;
    if (d == 0)
        return (complex_t){0, 0};
    return (complex_t){(a.r * b.r + a.i * b.i) / d, (a.i * b.r - a.r * b.i) / d};
}

/* --- FORMATTING (Compact ENG) --- */
void fmt_eng(float val, char *out)
{
    if (fabsf(val) < 1e-20)
    {
        strcpy(out, "0");
        return;
    }

    float log_v = log10f(fabsf(val));
    float exp_f = floorf(log_v);
    int eng_exp = (int)(floorf(exp_f / 3.0f) * 3.0f);
    float mantissa = val * powf(10.0f, (float)-eng_exp);

    if (fabsf(mantissa - roundf(mantissa)) < 1e-4)
    {
        sprintf(out, "%.0f", mantissa);
    }
    else
    {
        sprintf(out, "%.3f", mantissa);
    }

    if (eng_exp != 0)
    {
        char tmp[10];
        sprintf(tmp, "E%d", eng_exp);
        strcat(out, tmp);
    }
}

// Rectangular: A+Bi (Compacted spaces)
void fmt_rect(complex_t val, char *out)
{
    char r_str[20], i_str[20];
    int zr = (fabsf(val.r) < 1e-20);
    int zi = (fabsf(val.i) < 1e-20);

    if (zr && zi)
    {
        strcpy(out, "0");
        return;
    }

    if (zr)
    {
        fmt_eng(val.i, i_str);
        sprintf(out, "%si", i_str);
        return;
    }
    if (zi)
    {
        fmt_eng(val.r, r_str);
        sprintf(out, "%s", r_str);
        return;
    }

    fmt_eng(val.r, r_str);
    fmt_eng(fabsf(val.i), i_str);

    // Removed extra spaces to fit more on screen
    if (val.i < 0)
        sprintf(out, "%s-%si", r_str, i_str);
    else
        sprintf(out, "%s+%si", r_str, i_str);
}

// Polar (Rad)
void fmt_polar_rad(complex_t val, char *out)
{
    float r = sqrtf(val.r * val.r + val.i * val.i);
    float th = atan2f(val.i, val.r);
    char r_str[20];
    fmt_eng(r, r_str);
    sprintf(out, "%s<%.5fr", r_str, th); // Compacted " rad" to "r"
}

// Phasor (Deg)
void fmt_phasor_deg(complex_t val, char *out)
{
    float r = sqrtf(val.r * val.r + val.i * val.i);
    float th_rad = atan2f(val.i, val.r);
    float th_deg = th_rad * (180.0f / 3.14159265f);
    char r_str[20];
    fmt_eng(r, r_str);
    sprintf(out, "%s<%.4fd", r_str, th_deg); // Compacted " deg" to "d"
}

/* --- SCROLL PRINT HELPER --- */
void print_scrolled(const char *text, int row, int col, int offset)
{
    os_SetCursorPos(row, col);
    int len = strlen(text);
    if (offset >= len)
    {
        os_PutStrFull(" "); // End of string
    }
    else
    {
        os_PutStrFull(text + offset);
    }
}

/* --- PARSER --- */
const unsigned char *p_str;
float parse_expr();

void skip_spaces()
{
    while (*p_str == ' ')
        p_str++;
}

float parse_factor()
{
    skip_spaces();
    float val = 0;
    unsigned char token = *p_str;

    if (token == '(')
    {
        p_str++;
        val = parse_expr();
        if (*p_str == ')')
            p_str++;
        return val;
    }
    if (token == '-' || token == 0xB0)
    {
        p_str++;
        return -parse_factor();
    }

    // Constants & Funcs
    if (token == 0xAC)
    {
        p_str++;
        return 3.14159f;
    }
    if (token == 0xBB)
    {
        p_str++;
        return 2.71828f;
    }
    if (token == 0xBC)
    {
        p_str++;
        return sqrtf(parse_factor());
    }
    if (token == 0xC1)
    {
        p_str++;
        return sinf(parse_factor());
    }
    if (token == 0xC2)
    {
        p_str++;
        return cosf(parse_factor());
    }
    if (token == 0xC3)
    {
        p_str++;
        return tanf(parse_factor());
    }
    if (token == 0xC4)
    {
        p_str++;
        return logf(parse_factor());
    }
    if (token == 0xC5)
    {
        p_str++;
        return log10f(parse_factor());
    }

    int idx = 0;
    while ((token >= '0' && token <= '9') || token == '.' || token == 0x1B || token == 0xB0 || token == '-')
    {
        if (token == 0x1B)
            num_buf[idx++] = 'e';
        else if (token == 0xB0 || token == '-')
            num_buf[idx++] = '-';
        else
            num_buf[idx++] = (char)token;
        p_str++;
        token = *p_str;
        if (idx > 30)
            break;
    }
    num_buf[idx] = '\0';
    if (idx > 0)
        val = strtof(num_buf, NULL);
    return val;
}

float parse_pow()
{
    float val = parse_factor();
    skip_spaces();
    if (*p_str == '^' || *p_str == 0xF0)
    {
        p_str++;
        val = powf(val, parse_factor());
    }
    return val;
}

float parse_term()
{
    float val = parse_pow();
    skip_spaces();
    while (*p_str == '*' || *p_str == '/')
    {
        char op = *p_str++;
        float val2 = parse_pow();
        if (op == '*')
            val *= val2;
        else if (val2 != 0)
            val /= val2;
    }
    return val;
}

float parse_expr()
{
    float val = parse_term();
    skip_spaces();
    while (*p_str == '+' || *p_str == '-')
    {
        char op = *p_str++;
        float val2 = parse_term();
        if (op == '+')
            val += val2;
        else
            val -= val2;
    }
    return val;
}

float get_real_input(const char *prompt)
{
    memset(input_buf, 0, sizeof(input_buf));
    os_SetCursorPos(8, 0);
    os_PutStrFull("                ");
    os_SetCursorPos(9, 0);
    os_PutStrFull("                ");
    os_SetCursorPos(9, 0);
    os_PutStrFull(prompt);
    os_GetStringInput(input_buf, input_buf, 40);
    if (strlen(input_buf) == 0)
        return 0;
    p_str = (unsigned char *)input_buf;
    return parse_expr();
}

/* --- SOLVER (2D SCROLLABLE) --- */
void solve_and_display()
{
    os_ClrHome();
    os_PutStrFull("Solving...");

    // 1. Solve to Global X
    memcpy(M, A, sizeof(A));
    int i, j, k;
    for (i = 0; i < n; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            complex_t factor = c_div(M[j][i], M[i][i]);
            for (k = i; k <= n; k++)
            {
                M[j][k] = c_sub(M[j][k], c_mul(factor, M[i][k]));
            }
        }
    }
    for (i = n - 1; i >= 0; i--)
    {
        complex_t sum = {0, 0};
        for (j = i + 1; j < n; j++)
        {
            sum = c_add(sum, c_mul(M[i][j], X[j]));
        }
        X[i] = c_div(c_sub(M[i][n], sum), M[i][i]);
    }

    // 2. Interactive Scroll Loop
    int scroll_idx = 0; // Vertical scroll (which variables)
    int scroll_x = 0;   // Horizontal scroll (text offset)

    while (1)
    {
        os_ClrHome();
        os_FontSelect(os_SmallFont);

        int row = 0;
        for (int k = scroll_idx; k < scroll_idx + 2 && k < n; k++)
        {
            // Header X1:
            os_SetCursorPos(row, 0);
            sprintf(fmt_buf, "X%d:", k + 1);
            os_PutStrFull(fmt_buf);

            // Rect (SCROLLED)
            fmt_rect(X[k], fmt_buf);
            print_scrolled(fmt_buf, row, 4, scroll_x);

            // Polar (SCROLLED)
            fmt_polar_rad(X[k], fmt_buf);
            print_scrolled(fmt_buf, row + 1, 1, scroll_x);

            // Phasor (SCROLLED)
            fmt_phasor_deg(X[k], fmt_buf);
            print_scrolled(fmt_buf, row + 2, 1, scroll_x);

            // Separator
            os_SetCursorPos(row + 3, 0);
            os_PutStrFull("--------------------------");

            row += 4;
        }

        // Footer Instructions
        os_SetCursorPos(9, 0);
        os_PutStrFull("Arrows:Scroll Clr:Back");

        uint8_t key = 0;
        while ((key = os_GetCSC()) == 0)
            ;

        // Vertical Scroll
        if (key == sk_Down)
        {
            if (scroll_idx < n - 1)
                scroll_idx++;
        }
        if (key == sk_Up)
        {
            if (scroll_idx > 0)
                scroll_idx--;
        }

        // Horizontal Scroll
        if (key == sk_Right)
        {
            scroll_x++;
        }
        if (key == sk_Left)
        {
            if (scroll_x > 0)
                scroll_x--;
        }

        if (key == sk_Clear)
            return;
    }
}

/* --- MAIN --- */
int main(void)
{
    memset(A, 0, sizeof(A));
    os_ClrHome();
    os_FontSelect(os_SmallFont);

    float nv = get_real_input("Unknowns (2-5): ");
    n = (int)nv;
    if (n < 2)
        n = 2;
    if (n > 5)
        n = 5;

    while (1)
    {
        os_ClrHome();
        sprintf(fmt_buf, "--- SIZE: %d ---", n);
        os_SetCursorPos(0, 0);
        os_PutStrFull(fmt_buf);
        sprintf(fmt_buf, "<  EQUATION %d  >", cur_r + 1);
        os_SetCursorPos(1, 0);
        os_PutStrFull(fmt_buf);

        int j;
        for (j = 0; j <= n; j++)
        {
            os_SetCursorPos(j + 2, 0);
            if (j == cur_c)
                os_PutStrFull("> ");
            else
                os_PutStrFull("  ");

            if (j < n)
                sprintf(fmt_buf, "X%d: ", j + 1);
            else
                sprintf(fmt_buf, "Con: ");
            os_PutStrFull(fmt_buf);

            fmt_rect(A[cur_r][j], fmt_buf);
            os_PutStrFull(fmt_buf);
        }

        os_SetCursorPos(9, 0);
        os_PutStrFull("Enter=Edit Zoom=Solve");

        uint8_t key = 0;
        while ((key = os_GetCSC()) == 0)
            ;

        if (key == sk_Left)
            cur_r = (cur_r > 0) ? cur_r - 1 : n - 1;
        if (key == sk_Right)
            cur_r = (cur_r < n - 1) ? cur_r + 1 : 0;
        if (key == sk_Up)
            cur_c = (cur_c > 0) ? cur_c - 1 : n;
        if (key == sk_Down)
            cur_c = (cur_c < n) ? cur_c + 1 : 0;

        if (key == sk_Enter)
        {
            float r = get_real_input("Real: ");
            float i = get_real_input("Imag: ");
            A[cur_r][cur_c] = (complex_t){r, i};

            // Auto Advance
            cur_c++;
            if (cur_c > n)
            {
                cur_c = 0;
                cur_r++;
                if (cur_r >= n)
                    cur_r = 0;
            }
        }

        if (key == sk_Zoom || key == sk_Graph)
            solve_and_display();
        if (key == sk_Clear)
            break;
    }
    return 0;
}