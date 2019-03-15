import numpy as np
import os
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import matplotlib
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MaxNLocator
from matplotlib import rcParams
from collections import OrderedDict
import sys
import pdb


# basic set up
def plot_setting(fontsize=16):
    plt.rcParams['font.family']         = "Times New Roman"
    plt.rcParams['font.weight']         = 'bold'
    plt.rcParams['axes.labelsize']      = fontsize
    plt.rcParams['axes.labelweight']    = 'bold'
    plt.rcParams['xtick.labelsize']     = fontsize
    plt.rcParams['ytick.labelsize']     = fontsize
    matplotlib.rcParams.update({'font.size': fontsize,'legend.fontsize': fontsize-4})
    plt.close("all")
    return

def plot_legend():
    handles, labels         = plt.gca().get_legend_handles_labels()
    by_label                = OrderedDict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())

def pic_save(fname='XXX', outdir='.'):
    # plot_legend()
    plt.savefig(outdir+'/'+fname+'.pdf',dpi='figure',format='pdf')
    plt.close('all')


class FormatScalarFormatter(matplotlib.ticker.ScalarFormatter):
    def __init__(self, fformat="%1.1f", offset=True, mathText=True):
        self.fformat = fformat
        matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,
                                                        useMathText=mathText)
    def _set_format(self, vmin, vmax):
        self.format = self.fformat
        if self._useMathText:
            self.format = '$%s$' % matplotlib.ticker._mathdefault(self.format)


# plot functions

class plot_wavefrom_dataset(object):
    """plot the rfhv measurements and rf from real data.
    """
    # Plot for RFHV kernel
    def pic_frame(self, title="erzsol3"):
        plot_setting(fontsize=11)
        # w: 12 h:6
        fig                 = plt.figure(figsize=(12, 5))
        plt.suptitle(title, fontweight='bold')
        gs                  = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1], figure=fig)
        ax1                 = fig.add_subplot(gs[0,:])
        ax2                 = fig.add_subplot(gs[1,:])
        ax3                 = fig.add_subplot(gs[2,:])
        gs.tight_layout(fig, rect=[0.05, 0.05, 0.95, 0.95])
        return ax1, ax2, ax3, fig

    def plot_RF(self, rf, time, error=[], ax=[], xmax=10, ymax=0.6, ymin=-0.2, fname='syn_RF', outdir='.'):
        """
        plot R Receiver functions to explain
        ====================================================================================
        ::: input parameters :::
        time            - time (s) for the synthetic waveforms
        shift           - shift time for in deconvolution of RF (s)
        rf              - R component RF
        rfz             - Z component RF
        """
        if not ax:
            plt.close('all')
            self.plot_setting()
            fig                     = plt.figure(figsize=(7,5))
            ax                      = plt.subplot()
        # time                        = np.arange(rf.size)*dt
        # plot waveform
        ax.plot(time, rf,'r-', linewidth=1.5, label='true')
        if len(error):
            ax.fill_between(time, y2=rf-error, y1=rf+error, color='gray', alpha=0.5, lw=0, interpolate=True, label='error')
        # ax.plot([0,0],[0,max(rf)],'r--',lw=0.6)

        ax.set_xlim([0, xmax])
        ax.set_ylim([ymin, ymax])
        ax.xaxis.set_ticks(np.arange(0, xmax+2., step=2.))
        ax.yaxis.set_ticks(np.arange(ymin, ymax+0.2, step=0.2))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        # ax.set_xlabel('t (s)')
        ax.set_ylabel('amp', labelpad=-2)

        if not ax:
            self.plot_legend()
            plt.tight_layout()
            plt.savefig(outdir+'/'+fname+'.pdf',dpi='figure',format='pdf')
            plt.close('all')
        return

    def plot_wavefrom(self, waveform, time, ax=[], xmax=[], ymax=[], ymin=[], label='', fname='syn_RF', outdir='.'):
        """
        plot R Receiver functions to explain
        ====================================================================================
        ::: input parameters :::
        time            - time (s) for the synthetic waveforms
        shift           - shift time for in deconvolution of RF (s)
        rf              - R component RF
        rfz             - Z component RF
        """
        if not ax:
            plt.close('all')
            self.plot_setting()
            fig                     = plt.figure(figsize=(7,5))
            ax                      = plt.subplot()
        # time                        = np.arange(rf.size)*dt
        # plot waveform
        ax.plot(time, waveform,'k-', linewidth=1.5, label=label)

        if xmax:
            ax.set_xlim([0, xmax])
            ax.xaxis.set_ticks(np.arange(0, xmax+2., step=2.))

        # if ymax:
        #     ax.set_ylim([ymin, ymax])
        #     ax.yaxis.set_ticks(np.arange(ymin, ymax+0.2, step=0.2))

        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        # ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.set_xlabel('t (s)')
        ax.set_ylabel('amp', labelpad=2)
        ax.legend()
        if not ax:
            self.plot_legend()
            plt.tight_layout()
            plt.savefig(outdir+'/'+fname+'.pdf',dpi='figure',format='pdf')
            plt.close('all')
        return
