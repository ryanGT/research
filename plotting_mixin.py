class item_that_plots(object):
    def _get_fig(self, fi=1, fig=None):
            if fig is None:
                import pylab
                fig = pylab.figure(fi)
            return fig


    def _prep_ax(self, fi=1, fig=None, clear=True):
        fig = self._get_fig(fi=fi, fig=fig)
        if clear:
            fig.clf()
        self.ax = fig.add_subplot(1,1,1)
        return self.ax


    def label_axis(self, xlabel='Time (sec)', \
                   ylabel='Signal Amplitude (counts)'):
        self.ax.set_ylabel(ylabel)
        self.ax.set_xlabel(xlabel)

