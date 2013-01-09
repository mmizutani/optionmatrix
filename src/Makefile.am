bin_PROGRAMS = optionmatrix optionmatrix_console
noinst_PROGRAMS = checks

TESTS = checks

checks_SOURCES = \
checks.cpp\
args.cpp\
cumulative_norm.cpp\
defaults.cpp\
defs.h\
distributions.cpp\
extern.h\
futures.cpp\
greeks.cpp\
license.h\
logger.cpp\
main.h\
options_calls.cpp\
options_puts.cpp\
prototypes.h\
sanity.cpp\
termstructure.cpp\
time.cpp\
abradford/basicmodels.h\
abradford/black_76.cpp\
abradford/black_scholes.cpp\
abradford/delta.cpp\
abradford/future.cpp\
abradford/merton_73.cpp\
abradford/pre-BSM.cpp\
spinsky/binomial.cpp\
spinsky/binomial.h

optionmatrix_console_SOURCES = \
curses_calendar.cpp\
curses_demo.cpp\
curses_display.cpp\
curses_include.h\
curses_inputs.cpp\
curses_main.cpp\
curses_process.cpp\
curses_prop.cpp\
curses_misc.cpp\
args.cpp\
cumulative_norm.cpp\
defaults.cpp\
defs.h\
distributions.cpp\
extern.h\
futures.cpp\
greeks.cpp\
license.h\
logger.cpp\
main.h\
options_calls.cpp\
options_puts.cpp\
prototypes.h\
sanity.cpp\
termstructure.cpp\
tickers.cpp\
time.cpp\
abradford/basicmodels.h\
abradford/black_76.cpp\
abradford/black_scholes.cpp\
abradford/delta.cpp\
abradford/future.cpp\
abradford/merton_73.cpp\
abradford/pre-BSM.cpp\
spinsky/binomial.cpp\
spinsky/binomial.h

optionmatrix_SOURCES = \
gtk_bonds.cpp\
gtk_calendar.cpp\
gtk_cashflow.cpp\
gtk_checkbutton.cpp\
gtk_combobox.cpp\
gtk_dialog.cpp\
gtk_display.cpp\
gtk_extern.h\
gtk_glade.h\
gtk_include.h\
gtk_main.cpp\
gtk_main.h\
gtk_menu.cpp\
gtk_radio.cpp\
gtk_scale.cpp\
gtk_spin.cpp\
gtk_thread.cpp\
args.cpp\
cumulative_norm.cpp\
defaults.cpp\
defs.h\
distributions.cpp\
extern.h\
futures.cpp\
greeks.cpp\
license.h\
logger.cpp\
main.h\
options_calls.cpp\
options_puts.cpp\
prototypes.h\
sanity.cpp\
termstructure.cpp\
tickers.cpp\
time.cpp\
abradford/basicmodels.h\
abradford/black_76.cpp\
abradford/black_scholes.cpp\
abradford/delta.cpp\
abradford/future.cpp\
abradford/merton_73.cpp\
abradford/pre-BSM.cpp\
spinsky/binomial.cpp\
spinsky/binomial.h

optionmatrix_console_CPPFLAGS = -I$(top_srcdir)/src/financialrecipes -I$(top_srcdir)/src/metaoptions/src
optionmatrix_CPPFLAGS = -I$(top_srcdir)/src/financialrecipes -I$(top_srcdir)/src/metaoptions/src

optionmatrix_console_LDADD = ../src/financialrecipes/src/librecipes.a ../src/metaoptions/src/libmetaoptions.a
optionmatrix_LDADD = ../src/financialrecipes/src/librecipes.a ../src/metaoptions/src/libmetaoptions.a

checks_CPPFLAGS = -I$(top_srcdir)/src/financialrecipes -I$(top_srcdir)/src/metaoptions/src
checks_LDADD = ../src/financialrecipes/src/librecipes.a ../src/metaoptions/src/libmetaoptions.a

CLEANFILES = optionmatrix.log
