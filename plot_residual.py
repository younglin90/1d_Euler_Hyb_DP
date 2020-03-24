import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import style

#style.use('fivethirtyeight')


def animate(i):
    graph_data = open('residual.txt','r').read()
    lines = graph_data.split('\n')      # n 는 '엔터'의 뜻.
                                       # 파일 내의 '행'을 읽어서 lines에 저장.
    xs = []
    ys1 = []
    ys2 = []
    ys3 = []
    ys4 = []
    #ys5 = []
    for line in lines:
        if len(line) > 1:
            x,y1,y2,y3,y4 = line.split()   # split() 는 띄어쓰기로 구분된
            xs.append(float(x))           # 데이터들을 실수형식으로 추가해준다.
            ys1.append(float(y1))
            ys2.append(float(y2))
            ys3.append(float(y3))
            ys4.append(float(y4))
            #ys5.append(y5)

    ax1.clear()                 # 그래프들을 없애줌.
    ax2.clear()
    ax3.clear()
    ax4.clear()

    ax1.plot(xs,ys1,label='P',color='black',linewidth=0.5) # 플랏팅
    ax2.plot(xs,ys2,label='U',color='blue',linewidth=0.5)
    ax3.plot(xs,ys3,label='T',color='red',linewidth=0.5)
    ax4.plot(xs,ys4,label='R',color='orange',linewidth=0.5)

    ax1.set_title('P residual')
    ax2.set_title('U residual')
    ax3.set_title('R residual')
    ax4.set_title('T residual')

    ax1.set_xlabel('iteration')
    ax2.set_xlabel('iteration')
    ax3.set_xlabel('iteration')
    ax4.set_xlabel('iteration')

    ax1.set_ylabel('norm2')
    ax2.set_ylabel('norm2')
    ax3.set_ylabel('norm2')
    ax4.set_ylabel('norm2')

fig = plt.figure(figsize=(10,5))  # 플랏 토대 만들기
ax1 = fig.add_subplot(2,2,1)   # 서브플랏들 만들기
ax2 = fig.add_subplot(2,2,2)
ax3 = fig.add_subplot(2,2,3)
ax4 = fig.add_subplot(2,2,4)
fig.tight_layout(pad=4.0)    # 레이아웃들의 간격을 조절함.

ani = animation.FuncAnimation(fig,animate,interval=1000) # 실시간 그래프
plt.show()

