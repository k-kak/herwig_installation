BUILT_SOURCES = done-all-links

AUTOMAKE_OPTIONS = -Wno-portability

DIRLINKS = $(top_srcdir)/Amplitudes

CLEANFILES = done-all-links 

done-all-links: $(DIRLINKS)
	mkdir -p HJets
	$(LN_S) -f $(addprefix ../, $(DIRLINKS)) HJets
	touch done-all-links

install-data-local:
	find HJets -follow \( -name '*.h' -or -name '*.icc' \
	-or -name '*.tcc' -or -name '*.fh' -or -name '*.xh' \) \
          -exec $(install_sh_DATA) \{\} $(DESTDIR)$(includedir)/\{\} \;

uninstall-local:
	rm -rf $(DESTDIR)$(includedir)/HJets

clean-local:
	rm -rf HJets
