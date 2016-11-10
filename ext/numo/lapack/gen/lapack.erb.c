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

// Error Class ??
#define CHECK_DIM_GE(na,nd)                                     \
    if ((na)->ndim<(nd)) {                                      \
        rb_raise(nary_eShapeError,                              \
                 "n-dimension=%d, but >=%d is expected",        \
                 (na)->ndim, (nd));                             \
    }

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

<%
mod_var "mTL"

def_singleton("matmul",2,"matmul",mod_var:"mTL")
def_singleton("solve",2,"solve",mod_var:"mTL")
def_singleton("eigen",1,"eigen",mod_var:"mTL")

#def_alias "dot_mm", "matmul"

Function.codes.each do |x| %>
<%= x %>
<% end %>

void
Init_nary_<%=type_name%>_lapack()
{
    VALUE mN;
    mN = rb_const_get(rb_cObject, rb_intern("Numo"));
    cT = rb_const_get(mN, rb_intern("<%=class_name%>"));
    cCT = rb_const_get(mN, rb_intern("<%=complex_class_name%>"));
    mTL = rb_define_module_under(cT, "LAPACK");

    <% Function.definitions.each do |x| %><%= x %>
    <% end %>

    <% IdVar.assignment.each do |x| %>
    <%= x %><% end %>
}
