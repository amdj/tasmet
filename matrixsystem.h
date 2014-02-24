#ifndef MATRIXSYSTEM_H
#define MATRIXSYSTEM_H


class matrixsystem
{
	public:
		matrixsystem();

		dmat Kglob;
		dmat Rglob;

		us Ns;
		us TotalDofs;

		virtual ~matrixsystem();
	protected:

	private:
};

#endif // MATRIXSYSTEM_H
