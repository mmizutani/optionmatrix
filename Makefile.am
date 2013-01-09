SUBDIRS = src/financialrecipes/src src/financialrecipes/examples src/metaoptions/src src/metaoptions/tests src man doc doc/html

appicondir = $(datadir)/pixmaps
appicon_DATA = icons/32x32/optionmatrix.xpm\
	icons/32x32/optionmatrix_console.xpm

desktopdir = $(datadir)/applications
desktop_DATA = icons/optionmatrix.desktop\
	icons/optionmatrix_console.desktop

icondir = $(datadir)/icons/hicolor/256x256/apps
icon_DATA = icons/256x256/optionmatrix.png\
	icons/256x256/optionmatrix_console.png

icon48dir = $(datadir)/icons/hicolor/48x48/apps
icon48_DATA = icons/48x48/optionmatrix.png\
	icons/48x48/optionmatrix_console.png

icon32dir = $(datadir)/icons/hicolor/32x32/apps
icon32_DATA = icons/32x32/optionmatrix.png\
	icons/32x32/optionmatrix_console.png

EXTRA_DIST = \
	icons/32x32/optionmatrix.xpm\
	icons/32x32/optionmatrix_console.xpm\
	icons/32x32/optionmatrix.png\
	icons/32x32/optionmatrix_console.png\
	icons/48x48/optionmatrix.png\
	icons/48x48/optionmatrix_console.png\
	icons/256x256/optionmatrix.png\
	icons/256x256/optionmatrix_console.png\
	icons/optionmatrix.desktop\
	icons/optionmatrix_console.desktop\
	HOWTO\
	WEBSITE\
	src/gtkoptionmatrix.glade\
	src/images/10.png\
	src/images/11.png\
	src/images/1.png\
	src/images/2.png\
	src/images/3.png\
	src/images/4.png\
	src/images/5.png\
	src/images/6.png\
	src/images/7.png\
	src/images/8.png\
	src/images/9.png\
	src/images/installer.ico\
	src/images/optionmatrixdos.ico\
	src/images/optionmatrix.ico\
	utils/tab.sh\
	utils/push_to_web.sh\
	src/metaoptions/tests/TESTSUITE.README\
	src/metaoptions/AUTHORS\
	src/metaoptions/ChangeLog\
	src/financialrecipes/AUTHORS\
	src/financialrecipes/ChangeLog\
	src/financialrecipes/README\
	.gitignore

install: install-html install-pdf install-data
	$(MKDIR_P) $(DESTDIR)$(bindir)
	install src/optionmatrix_console $(DESTDIR)$(bindir)
	install src/optionmatrix $(DESTDIR)$(bindir)

website:
	sh utils/push_to_web.sh $(PACKAGE_URL) $(distdir).tar.gz $(PACKAGE).pdf
