/*
  lapack_<%=type_name%>.c
  Numerical Array Extension for Ruby
    (C) Copyright 1999-2007,2013 by Masahiro TANAKA

  This program is free software.
  You can distribute/modify this program
  under the same terms as Ruby itself.
  NO WARRANTY.
*/

#include <assert.h>
#include <ruby.h>
#include "numo/narray.h"
#include "numo/template.h"
//#include <blas.h>

#if defined __clang__
#  define UNUSED(name) __unused name
#else
#  define UNUSED(name) name
#endif

#define FFUNC(f) f##_

typedef <%=ctype%> dtype;
typedef <%=real_ctype%> rtype;
typedef <%=complex_type%> ctype;

<%
if is_complex
  def_id "real"
  def_id "imag"
end
%>
<%= IdVar.declaration.join("\n") %>

<% if is_complex %>
#include "numo/types/complex_macro.h"
<% end %>

static VALUE mTL;
static VALUE cT;
static VALUE cCT;
static VALUE cRT;
static VALUE cIT;

#define COUNT_OF_(a) (sizeof(a)/sizeof((a)[0]))

#define SET_POS(pos, i, type, n) do {(pos)[i] = (pos)[(i)-1] + ((sizeof(type)*(n)-1)/16+1)*16;} while (0)

// Error Class ??
#define CHECK_DIM_GE(na,nd) do {                                \
    if ((na)->ndim<(nd)) {                                      \
        rb_raise(nary_eShapeError,                              \
                 "n-dimension=%d, but >=%d is expected",        \
                 (na)->ndim, (nd));                             \
    }                                                           \
} while (0)

static size_t
max_(size_t m, size_t n)
{
  return ( (m > n) ? m : n ) ;
}

static size_t
min_(size_t m, size_t n)
{
  return ( (m < n) ? m : n ) ;
}

# if 0

// transpose copy
<% %w[_1d _2d].each do |d|
   %w[set get].each do |w| %>
static void
transpose_<%=w+d%>(dtype *x, na_loop_t *lp, int arg_idx, int n_args)
{
    size_t  i1, n1;
    ssize_t s1;
    char   *p1;
    <% if d=="_2d" %>
    size_t  i2, n2;
    ssize_t s2;
    char   *p2;
    <% end %>

    n1 = lp->args[arg_idx].shape[0]; // n[0+4*0] = n
    s1 = lp->args[arg_idx].iter[0].step;
    p1 = lp->args[arg_idx].ptr + lp->args[arg_idx].iter[0].pos;
    <% if d=="_2d" %>
    n2 = lp->args[arg_idx].shape[1]; // n[0+4*1] = n
    s2 = lp->args[arg_idx].iter[1].step;
    p2 = p1;
    for (i2=0; i2<n2; i2++) {
        p1 = p2;
        p2 += s2;
    <% else %>
    {
    <% end %>
        for (i1=0; i1<n1; i1++) {
            <% if w=="set" %>
            *(x++) = *(dtype*)(p1);
            <% else %>
            *(dtype*)(p1) = *(x++);
            <% end %>
            p1 += s1;
        }
    }
}
<% end end %>

#endif

<%
mod_var "mTL"

def_singleton("gemm",   2, "gemm",  mod_var:"mTL")
def_singleton("gesv",   2, "gesv",  mod_var:"mTL")
def_singleton("gels",   2, "gels",  mod_var:"mTL")
def_singleton("geqrf",  1, "geqrf", mod_var:"mTL")

def_singleton("getrf",  1, "getrf", mod_var:"mTL")
def_singleton("potrf", -1, "potrf", mod_var:"mTL")

def_singleton("gesvd", -1, "gesvd", mod_var:"mTL")
def_singleton("gesdd", -1, "gesdd", mod_var:"mTL")

def_singleton("matrix_rank", -1, "matrix_rank", mod_var:"mTL")  # this must be located after ges[vd]d

def_singleton("geev",  -1, "geev",  mod_var:"mTL")
def_singleton("heev",  -1, "heev",  mod_var:"mTL")
def_singleton("heevd", -1, "heevd", mod_var:"mTL")

def_singleton("norm", -1, "norm", mod_var:"mTL")

#def_alias "dot_mm", "matmul"

Function.codes.each do |x| %>
<%= x %>
<% end %>

void
Init_nary_<%=type_name%>_lapack()
{
    VALUE mN;
#if INT_MAX == 0x7fff
    char const * const int_name = "Int16";
#elif INT_MAX == 0x7fffffff
    char const * const int_name = "Int32";
#elif INT_MAX == 0x7fffffffffffffff
    char const * const int_name = "Int64";
#else
    char const * const int_name = 0;
#endif
    assert(int_name);
    mN = rb_const_get(rb_cObject, rb_intern("Numo"));
    cT = rb_const_get(mN, rb_intern("<%=class_name%>"));
    cCT = rb_const_get(mN, rb_intern("<%=complex_class_name%>"));
    cRT = rb_const_get(mN, rb_intern("<%=real_class_name%>"));
    cIT = rb_const_get(mN, rb_intern(int_name));
    mTL = rb_define_module_under(cT, "LAPACK");
    {
        size_t shape[1] = { 0 };
        VALUE const na = rb_narray_new(numo_cDFloat, COUNT_OF_(shape), &(shape[0]));
        ID const id_kahan_sum = rb_intern("kahan_sum");
        if (RTEST(rb_respond_to(na, id_kahan_sum))) {
            id_kahan_sum_or_sum = id_kahan_sum;
        } else {
            id_kahan_sum_or_sum = rb_intern("sum");
        }
    }

    <% Function.definitions.each do |x| %><%= x %>
    <% end %>

    <% IdVar.assignment.each do |x| %>
    <%= x %><% end %>
}
