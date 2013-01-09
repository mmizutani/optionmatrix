info_TEXINFOS = optionmatrix.texi

optionmatrix_TEXINFOS = optionmatrix.texi\
	fdl.texi\
	version.texi\
	epsf.tex

dvi: # do nothing for make dvi

CLEANFILES = optionmatrix.dvi

TEXTFILES = images/ConsoleChain.txt\
	images/ConsoleSelect.txt\
	images/ConsoleFormat.txt\
	images/ConsoleWelcome.txt\
	images/SleepDelay.txt\
	images/Calendar.txt\
	images/CouponEditor.txt\
	images/DateEngine.txt\
	images/OptionMatrixCalendarDated.txt\
	images/OptionMatrixSpreads.txt\
	images/OptionMatrix.txt\
	images/Time.txt\
	images/PriceInterestStandard.txt\
	images/StrikesDrop.txt\
	images/Strikes.txt\
	images/StrikeDropDownAndAdjuster.txt\
	images/ModelDropDownOpen.txt\
	images/CNDDrop.txt\
	images/CNDDropDown.txt\
	images/Greeks.txt\
	images/UseDateEngine.txt\
	images/DisplayFormats.txt\
	images/MonthsToStrike.txt\
	images/MonthToStrikes.txt\
	images/MonthsAcross.txt\
	images/3Calls3Puts.txt\
	images/CycleDropOpen.txt\
	images/DateEngineProperties.txt\
	images/LegControls.txt\
	images/CalendarSpread.txt\
	images/SpreadLegs.txt\
	images/PrefsSettings.txt\
	images/ModelsList.txt\
	images/Precision.txt\
	images/DateEngineSelected.txt\
	images/CycleClosed.txt\
	images/RealTime.txt\
	images/SpreadCheckbox.txt\
	images/SpreadsCCPPCPPC.txt\
	images/ConsoleCategories.txt\
	images/ConsoleInputs.txt\
	images/OptionMatrixWin1.txt\
	images/Leg1.txt\
	images/Leg1Leg2.txt\
	images/ModelByAuthor.txt\
	images/CalendarButton.txt\
	images/CUSTInput.txt\
	images/MonthAdjuster.txt\
	images/Export.txt

EQN_IMAGES = images/ConsoleChain.png\
	images/ConsoleSelect.png\
	images/ConsoleFormat.png\
	images/ConsoleWelcome.png\
	images/SleepDelay.png\
	images/CouponEditor.png\
	images/DateEngine.png\
	images/OptionMatrixCalendarDated.png\
	images/Calendar.png\
	images/OptionMatrix.png\
	images/OptionMatrixSpreads.png\
	images/3Calls3Puts.png\
	images/CalendarSpread.png\
	images/CNDDropDown.png\
	images/CNDDrop.png\
	images/CycleDropOpen.png\
	images/DateEngine.png\
	images/DateEngineProperties.png\
	images/DisplayFormats.png\
	images/Greeks.png\
	images/LegControls.png\
	images/ModelDropDownOpen.png\
	images/ModelsList.png\
	images/MonthsAcross.png\
	images/MonthsToStrike.png\
	images/MonthToStrikes.png\
	images/Precision.png\
	images/PrefsSettings.png\
	images/PriceInterestStandard.png\
	images/SpreadLegs.png\
	images/StrikeDropDownAndAdjuster.png\
	images/StrikesDrop.png\
	images/Strikes.png\
	images/Time.png\
	images/UseDateEngine.png\
	images/DateEngineSelected.png\
	images/CycleClosed.png\
	images/RealTime.png\
	images/SpreadCheckbox.png\
	images/SpreadsCCPPCPPC.png\
	images/ConsoleCategories.png\
	images/ConsoleInputs.png\
	images/OptionMatrixWin1.png\
	images/download.png\
	images/Leg1.png\
	images/Leg1Leg2.png\
	images/ModelByAuthor.png\
	images/CalendarButton.png\
	images/CUSTInput.png\
	images/MonthAdjuster.png\
	images/Export.png\
	images/6.png

EXTRA_DIST = $(EQN_IMAGES) $(TEXTFILES) html/images/README images/favicon.ico 

clean:
	-rm -f html/images/*.png
	-rm -f html/images/*.ico
	-rm -f html/*.html
	-rm -f html/optionmatrix.zip optionmatrix.pdf optionmatrix.info

install-info:
	install-info --info-dir="$(DESTDIR)$(infodir)/optionmatrix" \
		"$(DESTDIR)$(infodir)/optionmatrix.info";

install-html: $(optionmatrix__TEXINFOS) $(EQN_IMAGES)
	$(MKDIR_P) $(DESTDIR)${htmldir}/optionmatrix.html/images
	cp $(abs_top_srcdir)/src/images/1.png $(DESTDIR)${htmldir}/optionmatrix.html/images
	cp $(abs_top_srcdir)/src/images/2.png $(DESTDIR)${htmldir}/optionmatrix.html/images
	cp $(abs_top_srcdir)/src/images/3.png $(DESTDIR)${htmldir}/optionmatrix.html/images
	cp $(abs_top_srcdir)/src/images/4.png $(DESTDIR)${htmldir}/optionmatrix.html/images
	cp $(abs_top_srcdir)/src/images/5.png $(DESTDIR)${htmldir}/optionmatrix.html/images
	cp $(abs_top_srcdir)/src/images/6.png $(DESTDIR)${htmldir}/optionmatrix.html/images
	cp $(abs_top_srcdir)/src/images/7.png $(DESTDIR)${htmldir}/optionmatrix.html/images
	cp $(abs_top_srcdir)/src/images/8.png $(DESTDIR)${htmldir}/optionmatrix.html/images
	cp $(abs_top_srcdir)/src/images/9.png $(DESTDIR)${htmldir}/optionmatrix.html/images
	cp $(abs_top_srcdir)/src/images/10.png $(DESTDIR)${htmldir}/optionmatrix.html/images
	cp $(abs_top_srcdir)/src/images/11.png $(DESTDIR)${htmldir}/optionmatrix.html/images
	$(MAKEINFO) $(AM_MAKEINFOFLAGS) $(MAKEINFOFLAGS) -I $(srcdir) \
		--html --number-sections -o $(DESTDIR)${htmldir}/optionmatrix.html optionmatrix.texi
	for i in $(EQN_IMAGES); do cp -f ${srcdir}/$$i $(DESTDIR)${htmldir}/optionmatrix.html/images ; done

html: $(optionmatrix__TEXINFOS) $(EQN_IMAGES)
	$(MAKEINFO) $(AM_MAKEINFOFLAGS) $(MAKEINFOFLAGS) -I $(srcdir) \
		--html --number-sections -o html $(DESTDIR)optionmatrix.texi
	for i in $(EQN_IMAGES); do cp -f ${srcdir}/$$i $(DESTDIR)html/images ; done
	cp $(abs_top_srcdir)/src/images/1.png $(DESTDIR)html/images
	cp $(abs_top_srcdir)/src/images/2.png $(DESTDIR)html/images
	cp $(abs_top_srcdir)/src/images/3.png $(DESTDIR)html/images
	cp $(abs_top_srcdir)/src/images/4.png $(DESTDIR)html/images
	cp $(abs_top_srcdir)/src/images/5.png $(DESTDIR)html/images
	cp $(abs_top_srcdir)/src/images/6.png $(DESTDIR)html/images
	cp $(abs_top_srcdir)/src/images/7.png $(DESTDIR)html/images
	cp $(abs_top_srcdir)/src/images/8.png $(DESTDIR)html/images
	cp $(abs_top_srcdir)/src/images/9.png $(DESTDIR)html/images
	cp $(abs_top_srcdir)/src/images/10.png $(DESTDIR)html/images
	cp $(abs_top_srcdir)/src/images/11.png $(DESTDIR)html/images
	cp images/favicon.ico $(DESTDIR)html

maintainer-clean-local:
	-rm -rf html/images/*
	-rm -rf html/*.html
