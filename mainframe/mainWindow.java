package mainframe;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.StandardChartTheme;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.time.Millisecond;
import org.jfree.data.time.TimeSeries;
import org.jfree.data.time.TimeSeriesCollection;
import org.jfree.data.xy.XYDataset;

import com.sun.jna.Pointer;
import com.sun.jna.ptr.IntByReference;

public class mainWindow extends JFrame implements ActionListener, Runnable {
	/**
	 * 
	 */
	private static final long serialVersionUID = -666405807473668315L;
	static JButton start = new JButton("start");
	JProgressBar energy;
	private TimeSeries series;
	static public double lastValue = 100.0;
	static boolean flag = false;
	static boolean rec = true;

	void buildConstraints(GridBagConstraints gbc, int gx, int gy, int gw,
			int gh, int wx, int wy) {
		gbc.gridx = gx;
		gbc.gridy = gy;
		gbc.gridwidth = gw;
		gbc.gridheight = gh;
		gbc.weightx = wx;
		gbc.weighty = wy;
	}

	// layout
	@SuppressWarnings("deprecation")
	public mainWindow() {
		// frame
		super("注意力评估系统");
		setSize(1366, 720);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		JPanel pane = new JPanel();
		GridBagLayout grid = new GridBagLayout();
		GridBagConstraints constrains = new GridBagConstraints();
		pane.setLayout(grid);

		// energy bar
		energy = new JProgressBar();
		CreateProcess(energy);
		buildConstraints(constrains, 0, 0, 1, 1, 50, 80);
		grid.setConstraints(energy, constrains);
		constrains.fill = GridBagConstraints.BOTH;
		pane.add(energy);

		// dynamic timeseries
		this.series = new TimeSeries("注意力集中度", Millisecond.class);
		TimeSeriesCollection dataset = new TimeSeriesCollection(this.series);
		ChartPanel chartPanel = new ChartPanel(createChart(dataset));
		chartPanel.setPreferredSize(new java.awt.Dimension(500, 270));
		buildConstraints(constrains, 1, 0, 1, 1, 50, 80);
		grid.setConstraints(chartPanel, constrains);
		constrains.fill = GridBagConstraints.BOTH;
		pane.add(chartPanel);

		// start button
		start.addActionListener(this);
		buildConstraints(constrains, 0, 1, 2, 1, 100, 10);
		constrains.fill = GridBagConstraints.NONE;
		constrains.anchor = GridBagConstraints.CENTER;
		grid.setConstraints(start, constrains);
		pane.add(start);

		pane.setBackground(Color.WHITE);
		setContentPane(pane);
		setVisible(true);
	}

	// 创建能量条
	private void CreateProcess(JProgressBar energy) {
		energy.setOrientation(JProgressBar.VERTICAL);
		energy.setMaximum(200);
		energy.setMinimum(0);
		energy.setValue(0);
		energy.setBackground(Color.white);
		energy.setPreferredSize(new Dimension(100, 300));
		energy.setStringPainted(true);
		energy.setString(String.valueOf(lastValue));
	}

	// 创建时序图
	private JFreeChart createChart(XYDataset dataset) {
		// 创建主题样式
		StandardChartTheme standardChartTheme = new StandardChartTheme("CN");
		// 设置标题字体
		standardChartTheme.setExtraLargeFont(new Font("隶书", Font.BOLD, 20));
		// 设置图例的字体
		standardChartTheme.setRegularFont(new Font("宋书", Font.PLAIN, 15));
		// 设置轴向的字体
		standardChartTheme.setLargeFont(new Font("宋书", Font.PLAIN, 15));
		// 应用主题样式
		ChartFactory.setChartTheme(standardChartTheme);

		JFreeChart result = ChartFactory.createTimeSeriesChart("集中度趋势",
				"系统当前时间", "动态数值变化", dataset, true, true, false);
		XYPlot plot = (XYPlot) result.getPlot();
		ValueAxis axis = plot.getDomainAxis();
		axis.setAutoRange(true);
		axis.setFixedAutoRange(60000.0);
		axis = plot.getRangeAxis();
		axis.setRange(0.0, 200.0);
		return result;
	}

	/**
	 * 动态运行
	 */
	public void run() {
		// Thread thisThread = Thread.currentThread();
		while (true) {
			if (flag && !rec) {
				// double factor = 0.90 + 0.2 * Math.random();
				// lastValue = lastValue+1;
				Millisecond now = new Millisecond();
				this.series.add(now, lastValue);
				this.energy.setValue((int) lastValue);
				if(lastValue<90)
				energy.setForeground(Color.GREEN);
				else if(lastValue<130)
					energy.setForeground(Color.yellow);
				else if(lastValue<170)
					energy.setForeground(Color.CYAN);
				else energy.setForeground(Color.RED);
				energy.setString(String.valueOf((int) lastValue));
				rec = true;
			}
			try {
				Thread.currentThread();
				Thread.sleep(8);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
	}

	// event
	public void actionPerformed(ActionEvent evt) {
		Object src = evt.getSource();
		if (src == start) {
			String name = start.getText();
			if (name.equals("start")) {
				start.setText("stop");
				flag = true;
			} else {
				start.setText("start");
				flag = false;
			}
		}
	}

	public static void main(String[] args) {
		// create window
		mainWindow frame = new mainWindow();// TODO 自动生成的方法存根
		Thread map = new Thread(frame);
		map.start();

		// Receive EEgdata
		int recived_count = 0;
		Pointer eEvent = Edk.INSTANCE.EE_EmoEngineEventCreate();
		Pointer eState = Edk.INSTANCE.EE_EmoStateCreate();
		IntByReference userID = null;
		IntByReference nSamplesTaken = null;
		short composerPort = 1726;
		int option = 1;
		int state = 0;
		float secs = 1;
		int ED_AF3 = 3;
		int ED_AF4 = 16;
		int ED_F3 = 5;
		int ED_F4 = 14;
		boolean readytocollect = false;
//		double[] AF3 = new double[128];
//		double[] AF4 = new double[128];
//		double[] F3 = new double[128];
//		double[] F4 = new double[128];
		Complex[] AF3_tem = new Complex[128];
		Complex[] AF4_tem = new Complex[128];
		Complex[] F3_tem = new Complex[128];
		Complex[] F4_tem = new Complex[128];
		double[] m_AF3 = new double[32];//32
		double[] m_AF4 = new double[32];
		double[] m_F3 = new double[32];
		double[] m_F4 = new double[32];

		userID = new IntByReference(0);
		nSamplesTaken = new IntByReference(0);
		switch (option) {
		case 1: {
			if (Edk.INSTANCE.EE_EngineConnect("Emotiv Systems-5") != EdkErrorCode.EDK_OK
					.ToInt()) {
				System.out.println("Emotiv Engine start up failed.");
				return;
			}
			break;
		}
		case 2: {
			System.out.println("Target IP of EmoComposer: [127.0.0.1] ");

			if (Edk.INSTANCE.EE_EngineRemoteConnect("127.0.0.1", composerPort,
					"Emotiv Systems-5") != EdkErrorCode.EDK_OK.ToInt()) {
				System.out
						.println("Cannot connect to EmoComposer on [127.0.0.1]");
				return;
			}
			System.out.println("Connected to EmoComposer on [127.0.0.1]");
			break;
		}
		default:
			System.out.println("Invalid option...");
			return;
		}

		Pointer hData = Edk.INSTANCE.EE_DataCreate();
		Edk.INSTANCE.EE_DataSetBufferSizeInSec(secs);
		//System.out.print("Buffer size in secs: ");
		System.out.println(secs);

	//	System.out.println("Start receiving EEG Data!");
		// int sympol = 1;
		while (true) {
			/*//TEST
			   if(lastValue>=400||lastValue<=0){ sympol = -1*sympol; }
			   lastValue +=sympol; rec = false;
			 */// TEST END

			state = Edk.INSTANCE.EE_EngineGetNextEvent(eEvent);

			// New event needs to be handled
			if (state == EdkErrorCode.EDK_OK.ToInt()) {
				int eventType = Edk.INSTANCE.EE_EmoEngineEventGetType(eEvent);
				Edk.INSTANCE.EE_EmoEngineEventGetUserId(eEvent, userID);

				// Log the EmoState if it has been updated
				if (eventType == Edk.EE_Event_t.EE_UserAdded.ToInt())
					if (userID != null) {
						//System.out.println("User added");
						JOptionPane.showMessageDialog(null, "系统检测到已连接的Emotiv Epoc设备", "系统提示", JOptionPane.PLAIN_MESSAGE);
						Edk.INSTANCE.EE_DataAcquisitionEnable(
								userID.getValue(), true);
						readytocollect = true;
					}
			} else if (state != EdkErrorCode.EDK_NO_EVENT.ToInt()) {
				System.out.println("Internal error in Emotiv Engine!");
				break;
			}

			// 判断开始接收的两个条件，分别为设备状态和主窗体的开始按钮
		   if (readytocollect && flag) {
				Edk.INSTANCE.EE_DataUpdateHandle(0, hData);

				Edk.INSTANCE.EE_DataGetNumberOfSample(hData, nSamplesTaken);

				if (nSamplesTaken != null) {
					if (nSamplesTaken.getValue() != 0) {

						// System.out.print("Updated: ");
						// System.out.println(nSamplesTaken.getValue());

						double[] data = new double[nSamplesTaken.getValue()];
						if (recived_count <= 127) {// 如果接收满128个数据就停止接收
							for (int sampleIdx = 0; sampleIdx < nSamplesTaken
									.getValue(); ++sampleIdx) {
								if(recived_count>127) break;

								Edk.INSTANCE.EE_DataGet(hData, ED_AF3, data,
										nSamplesTaken.getValue());
								// 将接收到的数据转为复数以后复制到数组
								AF3_tem[recived_count] = new Complex(
										data[sampleIdx], 0);

								Edk.INSTANCE.EE_DataGet(hData, ED_AF4, data,
										nSamplesTaken.getValue());
								AF4_tem[recived_count] = new Complex(
										data[sampleIdx], 0);

								Edk.INSTANCE.EE_DataGet(hData, ED_F3, data,
										nSamplesTaken.getValue());
								F3_tem[recived_count] = new Complex(
										data[sampleIdx], 0);

								Edk.INSTANCE.EE_DataGet(hData, ED_F4, data,
										nSamplesTaken.getValue());
								F4_tem[recived_count] = new Complex(
										data[sampleIdx], 0);

								recived_count++;

							}
						}

					   if (recived_count > 127) {
							recived_count = 0;

							// 对四个通道的数据分别做傅里叶变换
							AF3_tem = fft(AF3_tem);
							AF4_tem = fft(AF4_tem);
							F3_tem = fft(F3_tem);
							F4_tem = fft(F4_tem);

							double value = 0;
							// 对傅里叶变换后的数据取幅值  j=7  j<12  j++
							for (int j = 7; j < 12; j++) {
								m_AF3[j] = power(AF3_tem[j]);
								m_AF4[j] = power(AF4_tem[j]);
								m_F3[j] = power(F3_tem[j]);
								m_F4[j] = power(F4_tem[j]);
								
								value = value+m_AF3[j]+m_AF4[j] +m_F3[j]+m_F4[j];
								/*
								if (rec) {
									// mainWindow.lastValue=AF3[j];
									mainWindow.lastValue = m_AF3[j]-3000 ;
									// System.out.println("rec" );
									rec = false;
								}
								*/
							}
							//value=1000-value;

							/// 非负矩阵分解
							pMatrix a = new pMatrix(new double[][] { m_AF3,m_AF4,
									m_F3, m_F4 });
							pMatrix[] b = a.nmf(16, 100, 0);
							//*/
							
							if (rec) {
								
								// mainWindow.lastValue=AF3[j];   Math.sqrt(b[0].stats()[2]);
								mainWindow.lastValue =  200-value*5;//放大100倍
								// System.out.println("rec" +mainWindow.lastValue  );
								
								rec = false;
							}
							
							/*
							System.out.println("a: " + a); // 原矩阵
							System.out.println("b[0]: " + b[0]); // 特征矩阵
							//System.out.println("b[0]第一行的和: " +  b[0].ValueOfRows(1,8,12)); 
							System.out.println("b[1]: " + b[1]); // 权系数矩阵
							//System.out.println("平方和: " + Math.sqrt(b[0].stats()[2]));
							pMatrix p = b[0].mult(b[1]); // 矩阵乘法
							System.out.println("b[0].mult(b[1]): " + p);
							double dist = a.dist(p); // 欧几里得距离
							System.out.println("dist: " + dist);
							*/
						}

					}
				}
			}
		}

		Edk.INSTANCE.EE_EngineDisconnect();
		Edk.INSTANCE.EE_EmoStateFree(eState);
		Edk.INSTANCE.EE_EmoEngineEventFree(eEvent);
		System.out.println("Disconnected!");

	}

	public static double power(Complex data) {
		double re = data.re();
		double im = data.im();
		double pow = Math.sqrt((re * re + im * im) / (128 * 128));
		return pow;
	}

	// compute the FFT of x[], assuming its length is a power of 2
	public static Complex[] fft(Complex[] x) {
		int N = x.length;

		// base case
		if (N == 1)
			return new Complex[] { x[0] };

		// radix 2 Cooley-Tukey FFT
		if (N % 2 != 0) {
			throw new RuntimeException("N is not a power of 2");
		}

		// fft of even terms
		Complex[] even = new Complex[N / 2];
		for (int k = 0; k < N / 2; k++) {
			even[k] = x[2 * k];
		}
		Complex[] q = fft(even);

		// fft of odd terms
		Complex[] odd = even; // reuse the array
		for (int k = 0; k < N / 2; k++) {
			odd[k] = x[2 * k + 1];
		}
		Complex[] r = fft(odd);

		// combine
		Complex[] y = new Complex[N];
		for (int k = 0; k < N / 2; k++) {
			double kth = -2 * k * Math.PI / N;
			Complex wk = new Complex(Math.cos(kth), Math.sin(kth));
			y[k] = q[k].plus(wk.times(r[k]));
			y[k + N / 2] = q[k].minus(wk.times(r[k]));
		}
		return y;
	}

}
